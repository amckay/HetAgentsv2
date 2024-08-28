
push!(LOAD_PATH,"./_aux")
using Parameters, ModelUtils, Plots
using LinearAlgebra: rank, I

include("HA_fcns.jl")
include("HANK_setup.jl")


compute_ss!(ss,par,grid);
T = 200;
print("Starting fake news...")
J = fakenews(ss,T,par,grid);
println("done")





# Now using ModelUtils to formulate the full model
print("Setting up for transition solve...")
@endogenousvariables Y π R i A τ
@exogenousvariables η ε

# in the HA block, we only take Jacobians for the relevant prices and there X = [W R]
# so we define an map from the full model to the ha block, this maps from the indices of the full X to the indices of the 
# HA X. If a variable is not in the HA X, we omit it
hamap_endog = Dict(1=>1,3=>2,6=>3);
hamap_exog = Dict(1=>4)

# steady state of aggregates --- note that we have one thing called ss and one called steady state, they are not the same
steadystate= Dict("initial" => [1.0;0.;par.R̄;par.R̄-1;par.Ā;ss[:X][3]], "exog" => [par.η̄;0.0]);

# Create model environment
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=T);

function fagg(m::ModelEnv,X,E)
    
    @unpack Y, π, R, i, A, τ = contemp(X,m)
    @unpack i_l, A_l  = lag(X,m)
    @unpack π_p = lead(X,m)
    
    @unpack η, ε = exogenous(E,m)

    @unpack κ, β, ϕ, Ā, η̄, R̄, ρA = m.par
    
    return [-R .+ (1 .+i_l)./(1 .+ π);# Fisher equation
            -π .+ κ .* log.(Y ./ m.ss["initial"][1]) .+ β * π_p;  # Phillips curve
            -i .+ (R̄-1).+ ϕ .* π .+ ε; # monetary policy rule
            -A .+  R.*A_l - τ .* Y .+ η; # govt budget
            -A .+ (1-ρA) *Ā .+ ρA .* A_l .+ η .- η̄ .+ (R.-1).*A_l .- (R̄-1)*Ā] # Treasury rule for τ
             
    
end
function fha(m::ModelEnv,X,E)
    @unpack Y,R,τ,A= contemp(X,m)
    @unpack η = exogenous(E,m)

    savings = solvesim(ss,[Y R τ η],par,grid);
    return savings .- A
end


fagg_X = get_fX(fagg,m);
fagg_E = get_fE(fagg,m);

fha_X = zeros(T,T*m.nx);
for i = 1:m.nx
    if haskey(hamap_endog,i)
        fha_X[:,(i-1)*T+1:i*T] = J[:,:,1,hamap_endog[i]];  # we use J[:,:,1,:] because we are using the first aggregator in the HA block, could do something else if you like
    end
end
# J has the Jacobian of solvesim, but fha is solvesim - A
fha_X[:,varindex(:A,m)] = -I(T);

fha_E = zeros(T,2T);  # the Jacobian of the HA block wrt exogenous variables
for i = 1:m.nexog
    if haskey(hamap_exog,i)
        fha_E[:,(i-1)*T+1:i*T] = J[:,:,1,hamap_exog[i]];  # we use J[:,:,1,:] because we are using the first aggregator in the HA block, could do something else if you like
    end
end



# Put it all together
f(m::ModelEnv,X,E) = [fagg(m,X,E);fha(m,X,E)];
checksteadystate(m,f,tol=5e-4)
f_X = [fagg_X;fha_X];
f_E = [fagg_E;fha_E];

println("done")



print("Linear transition..")
# setup the shock
Xss,E = longsteadystate(m);
E[varindex(:η,m,:exog)[1:4]] .= par.η̄ +0.0025; # fiscal shock


# linear solution
Xlin = Xss - f_X\f_E * E;
println("done")




# # nonlinear solution
print("Non-linear transition...")
MAXIT = 100
Xnl = copy(Xlin);
for it = 1:MAXIT
    println(it)
    global Xnl
    res = f(m,Xnl,E);
    if any(abs.(res) .> 1e-6)
        Xnl -= f_X \ res;
    else
        break
    end
end
display( plot(Xnl,m) )
savefig("../Figures/HANK.png")
