
using Parameters
using Roots: find_zero
include("_aux/rouwenhorst.jl")


# parameters
@with_kw mutable struct Par
    α = 0.4
    β = 0.98
    γ = 2.0
    δ = 0.02
    ρ = 0.95
    Ā = 4*1.4
    η̄ = 0.05
    R̄ = 1.005;
    κ = 0.1;
    ϕ = 1.5;
    ρA = 0.97;
end
par = Par();

struct HAGrids
    a::Array{Float64,1}
    na::Int64
    e::Array{Float64,1}
    ne::Int64
    Πe::Array{Float64,2}
end  

function HAGrids()

        


    amin = 0.0  # borrowing constraint
    #exogenous transition matrix
    #has form cols = future state, rows = current  state
    # outer constructor for CustomGrid
    ρskill = 0.978; 
    σskill = sqrt(0.01138); # See Floden and Linde note below
    Ne = 11;
    logendow, exogTrans = rouwenhorst(Ne,ρskill,σskill)
    endow = collect(exp.(logendow));

    # Floden-Linde
    # Their estimates are for an annual process.
    #     rhoQ = 0.9136^(1/4) = 0.978
    #
    #     same cross-sectional variance
    #     0.0426/(1-0.9136^2) = 0.2577 = sig2Q/(1-0.978^2)
    #     sig2Q = sig2A/(1-rhoA**2)*(1-rhoQ**2)
    #     sig2Q = 0.01138


    #labor endowments
    Lbar = sum(endow .* stationarydistribution(exogTrans'))
    endow = endow / Lbar
    Lbar = 1.0

    #grid for savings
    agridmin, agridmax, agridsize = amin, 400, 201
    agrid = range(agridmin^(0.25), agridmax^(0.25), length= agridsize).^4

    return HAGrids(agrid,agridsize,endow,Ne,exogTrans)

end
grid = HAGrids()

# functions for marginal utility and its inverse
uPrime = (par,c) -> c.^(-par.γ);
uPrimeInv = (par,up) -> up.^(-1.0/par.γ);

#initial guess of consumption function and initialize steady state dictionary
R = par.R̄
Y = 1;
c0 = (R-1)/R * grid.a .+ Y*grid.e';
τ =  ((R-1)*par.Ā + par.η̄)/Y;
ss = Dict{Symbol,Any}(:c =>c0,:X=>[Y R τ par.η̄])  # convention is X is a row vector


function aggregator(g,c,D)
    return D' * reshape(g,:);
end



"""
    EGMStepBack

    One step of the endogenous grid method.  Given the marginal value of assets
    Va_p in t+1 and prices at date t Xt, we produce the marginal value of assets at t 
    and the saving and consumption policy rules at t.
"""
function EGMStepBack(Va_p,Xt,par,grid)
    @unpack β, γ = par
    Y,R,τ,η = Xt

    uc = β * Va_p * grid.Πe';  # marginal utility today given Va' on grid of b'  
    cimplied = uPrimeInv(par, uc);   # convert to consumption
    
    currentassets = (grid.a .+ cimplied .- (1 .- τ) .* grid.e'*Y .- η  )/R; # current assets from budget constraint
    
    g = zeros(  size(currentassets)); 
    
    for j = 1:grid.ne
        g[:,j],_ = interp(currentassets[:,j],grid.a,grid.a,right_extrap=true); # interpolate savings rule onto grid
    end
    
    g = max.(g,grid.a[1])  # impose the borrowing constraint
    c = (R*grid.a .+ (1 .- τ) .* Y*grid.e' .+ η .-  g); # consumption from the budget constraint
    Va = R*uPrime(par,c); # envelope condition
    return Va, g, c
end


function checkβ!(β,par,grid) # this overwrites the global c0 and par.β
    # we fix the steady state real interest rate and calibrate β to clear the asset market
    par.β = β
    _, g, c0[:] = SolveEGM(c0,ss[:X],par,grid);
    Dss = stationarydistribution(maketransmat(g,grid));
    A_implied = aggregator(g,c0,Dss); 
    println("β = $β: market clearing = $(A_implied-par.Ā)")
    return A_implied - par.Ā
end


function compute_ss!(ss,par,grid)
    println("Computing steady state...")
    β_bracket = [0.97, 0.99]
    par.β = find_zero(b -> checkβ!(b,par,grid), β_bracket,atol=1e-7)
    ss[:Va], ss[:g], ss[:c] = SolveEGM(ss[:c],ss[:X],par,grid);  #compute the ss decision rules 
    ss[:D] = stationarydistribution(maketransmat(ss[:g],grid)); #compute the ss distribution
    ss[:agg] = aggregator(ss[:g],ss[:c],ss[:D]); # store the ss aggregates
    ss[:Nagg] = length(ss[:agg])  # it is convenient to have the size of the aggregates
    println("done")
end