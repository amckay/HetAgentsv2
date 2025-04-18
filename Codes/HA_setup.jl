
using Parameters
using Roots: find_zero

# parameters
@with_kw struct Par 
    α :: Float64 = 0.4
    β :: Float64 = 0.98
    γ :: Float64 = 2.0
    δ :: Float64 = 0.02
    ρ :: Float64 = 0.95
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
    Ne = 2
    lamw = 0.6 # job finding prob
    sigma = 0.2 # job separation prob
    exogTrans = [1-lamw  lamw; sigma  1.0-sigma]

    #labor endowments
    endow = [1.0;2.5]
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

function getWR(par,K;Z=1)
    @unpack α, δ = par
    Lbar = 1.0;
    W = (1-α)*Z.*(K /Lbar).^α;
    R = α*Z.*(K/Lbar).^(α-1) .+ 1 .- δ;
    return [W R]
end
#initial guess of consumption function and initialize steady state dictionary
KCompleteMarkets =  find_zero(K -> par.β * getWR(par,K)[2] - 1, 48.0) ;
W,R = getWR(par,48);
c0 = (R-1)/R * grid.a .+ W*grid.e';
ss = Dict{Symbol,Any}(:c =>c0)


function aggregator(g,c,D)
    return D' * reshape(g,:);
end


"""
    EGMStepBack

    One step of the endogenous grid method.  Given the marginal value of assets
    Va_p in t+1 and prices at date t Xt, we produce the marginal value of assets at t 
    and the saving and consumption policy rules at t.
"""
function EGMStepBack(Va_p::Array{Float64,2},Xt,par,grid)
    @unpack β, γ = par
    W,R = Xt

    uc = β * Va_p * grid.Πe';  # marginal utility today given Va' on grid of b'  
    cimplied = uPrimeInv(par, uc);   # convert to consumption
    
    currentassets = (grid.a .+ cimplied .- grid.e'*W)/R; # current assets from budget constraint
    
    g = zeros(  size(currentassets)); 
    
    for j = 1:grid.ne
        g[:,j],_ = interp(currentassets[:,j],grid.a,grid.a,right_extrap=true); # interpolate savings rule onto grid
    end
    
    g = max.(g,grid.a[1])  # impose the borrowing constraint
    c = (R*grid.a .+ W*grid.e' .-  g); # consumption from the budget constraint
    Va = R*uPrime(par,c); # envelope condition
    return Va, g, c
end



function checkK!(K,par,grid) # this overwrites the global c0
    _, g, c0[:] = SolveEGM(c0,getWR(par,K),par,grid);
    K_implied = aggregator(g,c0,stationarydistribution(maketransmat(g,grid))) 

    return K_implied - K
end


function compute_ss!(ss,par,grid)
    print("Computing steady state...")
    KLarge = 46.5;
    ss[:K] = find_zero(kk -> checkK!(kk,par,grid), (KCompleteMarkets+0.03, KLarge))
    ss[:X] = getWR(par,ss[:K]);  # store the steady state prices in :X
    ss[:Va], ss[:g], ss[:c] = SolveEGM(ss[:c],ss[:X],par,grid);  #compute the ss decision rules 
    ss[:D] = stationarydistribution(maketransmat(ss[:g],grid)); #compute the ss distribution
    ss[:agg] = aggregator(ss[:g],ss[:c],ss[:D]); # store the ss aggregates
    ss[:Nagg] = length(ss[:agg])  # it is convenient to have the size of the aggregates
    println("done")
end