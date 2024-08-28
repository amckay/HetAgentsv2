module ModelUtils

using ForwardDiff, Plots, SparseArrays
import Plots: plot
using LinearAlgebra: rank
using Symbolics
using Parameters

export VarList, ModelEnv, varindex, checksteadystate, longsteadystate, @endogenousvariables, @exogenousvariables, @addlaglead, contemp, lag, lead, exogenous, wide, long, get_fX, get_fE, linearIRFs, nonlineartransition, plot, optimaltransitionpath, optimalLQpolicy

abstract type VarList end


struct ModelEnv
    par
    vars::Dict{Any,DataType}
    T::Int
    ss::Dict{String,Vector{Float64}}
    nx::Int
    nexog::Int
end

function ModelEnv(;par,vars::Dict{Any,DataType},steadystate::Dict{String,Vector{Float64}},T::Int)


    nx = length(fieldnames(vars[0]));
    nexog = length(fieldnames(vars["exog"]));
    if ~haskey(steadystate,"terminal")
        steadystate["terminal"] = steadystate["initial"];
    end
    
    if length(steadystate["terminal"]) != nx || length(steadystate["initial"]) != nx
        error("Endogenous steady state is wrong size. Don't put commas in the variable list.")
    end
    if length(steadystate["exog"]) != nexog
        error("Exogenous steady state is wrong size. Don't put commas in the variable list.")
    end

    ModelEnv(par,vars,T,steadystate,nx,nexog)
end

"""return the position in X or E of variable s"""
function varindex(s::Symbol,m::ModelEnv,endogexog::Symbol=:endog)
    if endogexog == :endog
        fn = fieldnames(m.vars[0]);
    elseif endogexog == :exog
        fn = fieldnames(m.vars["exog"]);
    else
        error("argument not recognized")
    end
    j = findall(fn .== s);
    if length(j) == 1
        j = j[1];
    elseif length(j) == 0 
        error("$s not found")
    elseif length(j) > 1
        error("duplicate entries for $s")
    end
    return (j-1)*m.T+1:m.T*j
end


"""Check that the equations hold in steady state"""
function checksteadystate(m::ModelEnv,f;tol=1e-6)
    if any(abs.(m.ss["terminal"] .- m.ss["initial"]) .> tol)
        @warn "Cannot use checksteadystate when initial and terminal steady states are different"
        return
    end
    X, E = longsteadystate(m);
    res = f(m,X,E);
    if ~all(isapprox.(res, 0., atol=tol))
        println("Residuals in equations =$(findall(abs.( reshape(res,m.T,:)[1,:] )  .> tol))")
        error("Steady state does not hold")
    end

end

"""Macro that creates structs to reference endogenous variables and puts them in a dictionary called vardict"""
macro endogenousvariables(varnames...)

    
    e = quote 

        struct Vars <: VarList
            $([Symbol("$(entry)") for entry in varnames]...)
        end


        struct Vars_l <: VarList
            $([Symbol("$(entry)_l") for entry in varnames]...)
        end

        struct Vars_p <: VarList
            $([Symbol("$(entry)_p") for entry in varnames]...)
        end

        vardict = Dict{Any, DataType}(0=> Vars, -1=> Vars_l, 1=> Vars_p);
    end
    esc(e)

end

macro exogenousvariables(varnames...)
    e = quote
        
        struct VarsExog <: VarList
            $([Symbol("$(entry)") for entry in varnames]...)
        end

        vardict["exog"] = VarsExog;

    end

    esc(e)
end

macro addlaglead(t,varnames...)

    tt = eval(t)

    # build the  struct
    if tt > 0
        suff = "_p$(t)"
    else
        suff = "_l$(abs(tt))"
    end
    fields =[Symbol("$(entry)$(suff)") for entry in varnames]
    n2 = Symbol("Vars$(suff)")
    e = quote
        
        struct $(n2) <: VarList
        $(fields...)
        end

        vardict[$(t)] = $(n2);
    end
    esc(e)
end    


"""Return a VarList with the current period values"""
function contemp(X,m::ModelEnv)
    return m.vars[0](wide(X,m.nx)...)
end
"""Return a VarList with the lagged values"""
function lag(X,m::ModelEnv)
    return m.vars[-1](lagshift(wide(X,m.nx),m)...)
end
"""Return a VarList with the t-j values"""
function lag(X,m::ModelEnv,j::Int)
    return m.vars[-j](lagshift(wide(X,m.nx),m,j)...)
end
"""Return a VarList with the leaded values"""
function lead(X,m::ModelEnv)
    return m.vars[1](leadshift(wide(X,m.nx),m)...)
end
"""Return a VarList with the t+j values"""
function lead(X,m::ModelEnv,j::Int)
    return m.vars[j](leadshift(wide(X,m.nx),m,j)...)
end
"""Return a VarList with the exogenous values"""
function exogenous(E,m::ModelEnv)
    return m.vars["exog"](wide(E,m.nexog)...)
end

"""Take a vector x of vectors xi and lag each one by prepending the steady state value"""
function lagshift(x,m::ModelEnv)
    return [[m.ss["initial"][i]; xi[1:end-1]]  for (i,xi) in enumerate(x)]
end
"""Take a vector x of vectors xi and lag each j times one by prepending the steady state value"""
function lagshift(x,m::ModelEnv,j::Int)
    return [[repeat( [m.ss["initial"][i]] ,j); xi[1:end-j]]  for (i,xi) in enumerate(x)]
end
"""Take a vector x of vectors xi and lead each one by shifting and appending the steady state value"""
function leadshift(x,m::ModelEnv)
    return [[xi[2:end]; m.ss["terminal"][i]]  for (i,xi) in enumerate(x)]
end
"""Take a vector x of vectors xi and lead each one j times by shifting and appending the steady state value"""
function leadshift(x,m::ModelEnv,j::Int)
    return [[xi[j+1:end]; repeat( [m.ss["terminal"][i]] ,j)]  for (i,xi) in enumerate(x)]
end

"""Take an n * T vector and reshape it to n vectors of length T"""
function wide(x,n::Int)
    T =  Int(length(x)/n);
    return [x[1+T*(i-1):T*i] for i in 1:n]
end
"""Take n vectors of length T and reshape it to  n * T vector"""
function long(x)
    return vcat(x...)
end

"""Return a T*nx vector and T*nexog vector of steady state values"""
function longsteadystate(m::ModelEnv;ss = m.ss["terminal"])
    return reshape( ones(m.T) .* ss', m.T * m.nx),  reshape( ones(m.T) .* m.ss["exog"]', m.T * m.nexog)
end



"""Jacobian of f w.r.t. X evaluated at (X,E)
If (X,E) are not provided, we assume steady state
"""
function get_fX(f,m,X,E)
    return ForwardDiff.jacobian(XX -> f(m,XX,E),X)
end
function get_fX(f,m)
    return get_fX(f,m,longsteadystate(m)...)
end

"""Jacobian of f w.r.t. E evaluated at (X,E)

If (X,E) are not provided, we assume steady state
"""
function get_fE(f,m,X,E)
    return ForwardDiff.jacobian(EE -> f(m,X,EE),E)
end
function get_fE(f,m)
    return get_fE(f,m,longsteadystate(m)...)
end

"""
Get the linear IRFs akin to a first-order approximation

IRFMat = linearIRFs(f,m,X,E)

If X and E are provided, the approximation is taken around that path
X and E default to steady state.

IRFMat is a matrix nx*T x nexog*T giving the IRFs of each variable to each shock at each horizon
IRFMat is in deviations from X

"""    
function linearIRFs(f,m::ModelEnv,X,E;)
    fX  = get_fX(f,m,X,E);
    fE  = get_fE(f,m,X,E);
    return -sparse(fX)\fE # converting to sparse matrix means faster solve
end
function linearIRFs(f,m::ModelEnv)
    return linearIRFs(f,m,longsteadystate(m)...)
end



"""Solve for nonlinear transition path using X0 as initial guess"""
function nonlineartransition(m,f,X0,E;maxit=100)
    println("nonlinear transition solver:")
    X = copy(X0);
    for it = 1:maxit
        Xold = copy(X);
        res = f(m,X,E);
        fX = sparse(get_fX(f,m,X,E));
        X -= fX\res;
        update = maximum(abs.(X .- Xold));
        println("it = $it,  update = $update")
        if update < 1e-6
            return X
        end
    end
    error("non-linear solve did not converge")
end



"""Plot IRFs

IRFMat -- T*nx x T*nexog array with IRFs in the columns
m  -- model environment, which gives variable names
irf_horizon -- how many periods to plot
shock --  int in 1:nexog or symbol that indicates which shock to plot
shock_horizon -- by default we plot a date-0 shock (first column), but you could ask for an anticipated shcok
series_names-- by default we plot all the endogenous variables. But you could supply an array of symbols and then only those are plotted.
"""
function plot(IRFMat::Array{Float64,2},m::ModelEnv;irf_horizon=20,shock::Any=1,shock_horizon=0,series_names=fieldnames(m.vars[0]))
    
    if isa(shock,Int)
        shock_number = shock;
    elseif isa(shock,Symbol)
        shock_number = findall( fieldnames(m.vars["exog"]) .== shock);
        @assert(length(shock_number) == 1)
        shock_number = shock_number[1]
    else
        error("Invalid argument --- shock must be Int or Symbol")
    end


    IRFvec = IRFMat[:,(shock_number-1)*m.T+shock_horizon+1];
    return plot(IRFvec,m,irf_horizon=irf_horizon,series_names=series_names)    
end



"""Plot IRFs

IRFvec --  T*nx vector with IRFs to one shock
m  -- model environment, which gives variable names
irf_horizon -- how many periods to plot
series_names-- by default we plot all the endogenous variables. But you could supply an array of symbols and then only those are plotted.
"""
function plot(IRFvec::Vector{Float64},m::ModelEnv;irf_horizon=20,series_names=fieldnames(m.vars[0]))    

    IRFs = [contemp(IRFvec,m) ];
    plot(IRFs,m,irf_horizon=irf_horizon,series_names=series_names,shock_labels="")
end

"""
Plot a vector of IRFs stored in Vars objects

IRFs vector of Vars that each contain the IRFs to a particular shock
m  -- model environment, which gives variable names
irf_horizon -- how many periods to plot
series_names-- by default we plot all the endogenous variables. But you could supply an array of symbols and then only those are plotted.
shock_labels -- labels for the elements of IRFs that go in the legend
"""
function plot(IRFs::Vector{T} where T <: VarList,m::ModelEnv;irf_horizon=20,series_names=fieldnames(m.vars[0]),shock_labels="")    


    pltarr = Array{Plots.Plot,1}(undef,length(series_names));
    for (i,ser) in enumerate(series_names)

        dat = hcat([getfield(IRF,ser)[1:irf_horizon+1] for IRF in IRFs]...);
        
        if maximum(dat) - minimum(dat) < 1e-6
            ylim = dat[1] .+[-0.01;0.01];
        else
            ylim = :auto;
        end
        if i == 1 && size(dat,2) > 1
            label = shock_labels
        else
            label = ""
        end
        pltarr[i] = plot(0:irf_horizon, dat,label=label,ylims=ylim);
        title!(String(ser))
    end
    plot(pltarr...,size=(800,600))
end
    



# ----- optimal policy methods -----








"""Computes the Hessian of λ' * f"""
function _getlamfhessian(f,m,X,E,lam)
    @unpack nx,nexog,T = m
    nf = m.nx-1;

    Tsym = 5;
    @variables xs[1:nx*Tsym] lams[1:nf*Tsym] es[1:nexog*Tsym]
    H = Symbolics.hessian(sum(lams .* f(m,xs,es)),[xs...])

   
    HH = zeros(nx*T,nx*T);

    lamwide = reshape(lam,:,nf);
    Xwide = reshape(X,:,nx);
    Ewide = reshape(E,:,nexog);

    myconvert(n::Num) = Float64(n.val);
    for t = 2:T-1
        if t == 2
            ll = reshape([zeros(1,nf) ; lamwide[1:t+2,:]],:);
            xx = reshape([m.ss["initial"]' ; Xwide[1:t+2,:]],:);
            ee = reshape([ m.ss["exog"]'; Ewide[1:t+2,:]],:);
        elseif t == T-1
            ll = reshape([lamwide[t-2:T,:];lamwide[[T],:]],:);
            xx= reshape([Xwide[t-2:T,:]; m.ss["terminal"]'],:);
            ee = reshape([Ewide[t-2:T,:]; m.ss["exog"]'],:);
        else
            ll = reshape(lamwide[t-2:t+2,:],:);
            xx = reshape(Xwide[t-2:t+2,:],:);
            ee = reshape(Ewide[t-2:t+2,:],:);
        end
        
        Htemp = Symbolics.substitute(H,Dict(xs=>xx,lams=>ll,es=>ee));
        for i = 1:nx
            ofsi = (i-1)*T;
            for j = 1:nx
                ofsj = (j-1)*T;
                HH[ofsi .+ (t-1:t+1),ofsj .+ (t-1:t+1)] = myconvert.(Htemp[(i-1)*Tsym.+(2:4),(j-1)*Tsym.+(2:4)]);
            end
        end

    end
    return HH
end

"""optimaltransitionpath

Compute transition path with optimal policy.

Inputs: 
Ugrad  function of m,X,E that returns nX vector giving gradient of U
Uhess  function of m,X,E that gives nX x nX array giving the Hessian of U
f  function of m,X,E  giving private sector residuals at X,E
m model environment
X0 initial guess of transition path
E exogenous terms on the transition path

Returns:
X     optimal transition path
lam   Lagrange multipliers at the optimum
"""
function optimaltransitionpath(Ugrad,Uhess,f,m,X0,E; maxit=100, lam0=nothing)

    X = copy(X0);
    fX = get_fX(f,m,X,E);
    neq = size(fX,1);
    UX = Ugrad(m,X,E);

    #initial guess of lam 
    if isnothing(lam0)
        lam = fX' \ UX; 
    else
        lam = lam0;
    end

    for it = 1:maxit
        psres = f(m,X,E);
        res = [UX .- fX' * lam  ; psres];

        jac = [Uhess(m,X,E) .- _getlamfhessian(f,m,X,E,lam)     -fX';   fX   zeros(neq,neq) ];
        update = -jac \ res;

        X = X + update[1:m.nx*m.T];
        lam = lam + update[m.nx*m.T+1:end];

        fX = get_fX(f,m,X,E);
        UX = Ugrad(m,X,E);

        println("optimaltransitionpath: it = $it    update = $(maximum(abs.(update)))")
        if maximum(abs.(update)) < 1e-7
            break
        end

    end

    # if we are doing the non-linear transition, check residuals are zero
    if maxit > 1
        res = [UX .- fX' * lam  ; f(m,X,E)];
        if maximum(abs.(res[1:size(X,1)])) > 1e-6
            error("FOCs are not zero")
        end
        if maximum(abs.(res[size(X,1)+1:end])) > 1e-6
            error("Model equations are not zero")
        end
    end
    return X, lam
end

"""LQ policy
Optimal policy for an LQ problem.
  Objective: 0.5 X' Q X + X' R 
  Constraints: A + B X = 0
"""
function optimalLQpolicy(Q::AbstractArray,R::Array{Float64,1},A::Array{Float64,1},B::Array{Float64,2})
    nX = size(Q,1);
    @assert size(Q,2) == nX;
    @assert size(R) == (nX,)
    nf = size(B,1);
    @assert size(B,2) == nX
    @assert size(A) == (nf,)

    xlam = -[Q -B'; B zeros(nf,nf)] \ [R;A];

    return xlam[1:nX], xlam[nX+1:end]
end

"""LQ policy
Optimal policy for an LQ problem.
  Objective: U(X̂)≈Ugrad(X)(X̂ - X) + 0.5 (X̂ -X)' Uhess (X̂ -X)
  Constraints: f(X̂) ≈ f(X) + f_X(X) (X̂-X)  +  f_E Ê=  0
"""
function optimalLQpolicy(m::ModelEnv,Ugrad,Uhess,f,Ehat;approx_point=nothing)
    if isnothing(approx_point)
        X,E = longsteadystate(m);
    else
        X,E = approx_point;
    end
    R = Ugrad(m,X,E);
    Q = Uhess(m,X,E);
    B = get_fX(f,m,X,E);
    A = get_fE(f,m,X,E) * Ehat .+ f(m,X,E) ;

    dX, lam = optimalLQpolicy(Q,R,A,B);
    return X .+ dX, lam

end

end # end module