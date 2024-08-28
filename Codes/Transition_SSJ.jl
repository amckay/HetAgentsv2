
push!(LOAD_PATH,"./_aux")
using Parameters, ModelUtils, Plots
using LinearAlgebra: rank, I

include("HA_fcns.jl")
include("HA_setup.jl")


compute_ss!(ss,par,grid);
T = 200;
print("Starting fake news...")
J = fakenews(ss,T,par,grid)
println("done")





# Now using ModelUtils to formulate the full model
print("Setting up for transition solve...")
@endogenousvariables K W R
@exogenousvariables Z

# in the HA block, we only take Jacobians for the relevant prices and there X = [W R]
# so we define an map from the full model to the ha block, this maps from the indices of the full X to the indices of the 
# HA X. If a variable is not in the HA X, we omit it
hamap = Dict(2=>1,3=>2);

# steady state of aggregates --- note that we have one thing called ss and one called steady state, they are not the same
steadystate= Dict("initial" => [ss[:K]; ss[:X]...], "exog" => [1.]);

# Create model environment
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=T);

function fagg(m::ModelEnv,X,E)
    
    @unpack K, W, R = contemp(X,m)
    @unpack K_l = lag(X,m)
    
    @unpack Z = exogenous(E,m)

    @unpack α, β, γ, δ, ρ = m.par
    
    return [-R .+ α * Z .* K_l.^(α-1) .+ 1 .- δ;
            -W .+ (1-α) .*Z .*K_l.^α ]
    
end
function fha(m::ModelEnv,X)
    @unpack K, W, R = contemp(X,m)

    savings = solvesim(ss,[W R],par,grid);
    return savings .- K 
end


 # the residuals in the HA equations tend to be larger because finding the steady state involves many computations

fagg_X = get_fX(fagg,m);
fagg_E = get_fE(fagg,m);

fha_X = zeros(T,T*m.nx);
for i = 1:m.nx
    if haskey(hamap,i)
        fha_X[:,(i-1)*T+1:i*T] = J[:,:,1,hamap[i]];  # we use J[:,:,1,:] because we are using the first aggregator in the HA block, could do something else if you like
    end
end
# J has the Jacobian of solvesim, but fha is solvesim - K
fha_X[:,varindex(:K,m)] = -I(T);

fha_E = zeros(T,T);  # the exogenous variable does not direclty enter the HA block




# Put it all together
f(m::ModelEnv,X,E) = [fagg(m,X,E);fha(m,X)];
checksteadystate(m,f,tol=5e-4)
f_X = [fagg_X;fha_X];
f_E = [fagg_E;fha_E];

println("done")



print("Linear transition..")
# setup the shock
Xss,_ = longsteadystate(m);
Zpath = exp.( 0.03 * 0.9.^(0:T-1)  );


# linear solution
Xlin = Xss - f_X\f_E * (Zpath .- 1.0);
println("done")




# # nonlinear solution
print("Non-linear transition...")
MAXIT = 100
Xnl = copy(Xlin);
for it = 1:MAXIT
    println(it)
    global Xnl
    res = f(m,Xnl,Zpath);
    if any(abs.(res) .> 1e-6)
        Xnl -= f_X \ res;
    else
        break
    end
end
display( plot(Xnl,m) )
eqmK = Xnl[varindex(:K,m)];
println("done")
println("K for t = 0 to 4: $(eqmK[1:5])")