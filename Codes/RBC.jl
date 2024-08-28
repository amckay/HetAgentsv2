
push!(LOAD_PATH,"./_aux")
using Parameters, ModelUtils


# parameters
@with_kw struct Par 
    α = 0.4
    β = 0.98
    γ = 2.0
    δ = 0.02
    ρ = 0.95
end
par = Par();


@endogenousvariables Z R K Y C
@exogenousvariables ε



# steady state
function get_ss(par)
    @unpack α, β, γ, δ = par
    Z = 1.
    R = 1 / β
    K = ((R-1+δ)/α)^(1/(α-1))
    Y = K^α
    C = Y - δ*K
    return [Z; R; K; Y; C]
end
steadystate= Dict("initial" => get_ss(par), "exog" => [0.]);

# Create model environment
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=300);

function f(m::ModelEnv,X,E)
    
    @unpack Z, R, K, Y, C = contemp(X,m)
    @unpack C_p, R_p = lead(X,m)
    @unpack K_l, Z_l = lag(X,m)
    
    @unpack ε = exogenous(E,m)

    @unpack α, β, γ, δ, ρ = m.par
    
    return [-C.^(-γ) .+ β * R_p .* C_p.^(-γ);
            -R .+ α * Z .* K_l.^(α-1) .+ 1 .- δ;
            -K .+ (1-δ)*K_l .+ Y .- C;
            -Y .+ Z .* K_l.^α;
            -log.(Z) .+ ρ .* log.(Z_l) .+ ε]
    
end

checksteadystate(m,f)

# setup the shock
Xss,E = longsteadystate(m);
E[1] = 0.03; # 3% shock to TFP at date 0


# linear solution
Xlin = Xss + linearIRFs(f,m) * E;
display( plot(Xlin,m) )

# nonlinear solution
Xnl = nonlineartransition(m,f,Xlin,E);
display( plot(Xnl,m) )