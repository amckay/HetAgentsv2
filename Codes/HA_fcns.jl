using Parameters
using LinearAlgebra: eigen


function stationarydistribution(T;method=:eigen,maxit=100_000)
    # columns of T should sum to 1
    @assert all(abs.(sum(T,dims=1)[:].-1) .<1e-12)
    if method==:eigen
        f = eigen(T);
        i = argmin(abs.(f.values .- 1));
        D = real(f.vectors[:,i]);
        return D / sum(D);
    elseif method==:simulate 
        D = zeros(size(T,1))
        D[1] = 1.0;
        for it = 1:maxit
            D = T * D;
            if it % 100 == 0 &&  all( abs.(D - T*D) .< 1e-8)
                break
            end
        end
        return D
    else
        error("Method not recognized")
    end
end


"""
Linear Interpolation

    Arguments:
    x,y  vectors that define the function
    x1   vector of points to evaluate the function

    Returns:
    y1   vector of function values at x1

    When x1 is outside the range of x, we use the following approach
    if right_extrap == false, x1 > x[end] maps to y[end]
    if right_extrap == true, x1 > x[end] is extrapolated based on the last two points of y
    The treatment of x1< x[1] is analogous and depends on left_extrap.
    Both extrapolations default to false
"""    
function interp(x,y,x1;right_extrap=false, left_extrap=false)
    i = lookup(x,x1);
    xl = x[i.-1]
    xr = x[i]
    yl = y[i.-1]
    yr = y[i]
    y1 = yl .+ (yr .-yl)./(xr .-xl) .* (x1 .-xl)
    
    
    if ~right_extrap 
        above = x1 .> x[end];
        y1[above] .= y[end];
    end
    if ~left_extrap
        below = x1 .< x[1];
        y1[below] .= y1[1];
    end

    return y1, i
end

"""lookup
    x is a vector of sorted values
    x1 is a vector
    returns the indices if the right side of the bracket x_{i-1} < x <= x_{i}
"""
function lookup(x,x1)
    N = length(x)
    i = zeros(Int32,length(x1));
    for (j,x1j) in enumerate(x1)
        i[j] = min.(max.(searchsortedlast(x,x1j,lt= <=),1),N-1) .+1;
    end
    return i
end



"""SolveEGM

    Iterate the endogenous grid method to convergence
"""
function SolveEGM(c0,Xt,par,grid)
    Va = uPrime(par,c0);
    #loop until convergence
    # println("solving for policy rules")
    tol = 1e-10
    
    for it =1:10000
        Va1 = EGMStepBack(Va,Xt,par,grid)[1]

        if it % 50 == 0
            test = abs.(Va1-Va)/(abs.(Va) .+tol)
            # println("it = $(it), test = $(maximum(test))")
            if all(test  .< tol)
                break
            end
        end

        Va = Va1
    end
    return EGMStepBack(Va,Xt,par,grid)
end


function maketransmat(g,grid)
    gcapped = min.(g, grid.a[end]);
    # Rows of M correspond to where we are going, cols correspond to where we are coming from
    M = zeros(grid.ne*grid.na,grid.ne*grid.na)
    for j in 1:grid.ne
        i = lookup(grid.a,gcapped[:,j])
        p = (gcapped[:,j] .- grid.a[i.-1]) ./ (grid.a[i] .- grid.a[i.-1]);
        sj = (j-1)*grid.na
        for k =1:grid.ne
            sk = (k-1) * grid.na
            for m = 1:grid.na
                M[sk .+ i[m], sj+m] = p[m] * grid.Πe[j,k];
                M[sk .+ i[m] - 1,sj+m] = (1.0 .-p[m]) * grid.Πe[j,k];
            end
        end
    end
    @assert sum(M,dims=1)[:] ≈ ones(grid.ne*grid.na)
    
    return M
end

"""solvesim  -- For given path of prices (Xpath), solve backwards and simulate forwards
"""
function solvesim(ss,Xpath,par,grid)
    T = size(Xpath,1);
    gpols = zeros(size(ss[:g])...,T);
    cpols = zeros(size(ss[:c])...,T);

    nagg = size(aggregator(ss[:g],ss[:c],ss[:D]));# here we call aggregator at steady state just to see what size we need
    agg  = zeros(T,nagg...);  

    # Solve backwards in time to get the savings policy rules
    Va_p = copy(ss[:Va]);
    for t = T:-1:1
        Va_p,gpols[:,:,t],_ = EGMStepBack(Va_p,Xpath[t,:],par,grid);
    end


    # Simulate forward in time to get the distribution of wealth
    D = copy(ss[:D]);
    for t = 1:T
        agg[t,:] .= aggregator(gpols[:,:,t],cpols[:,:,t],D);
        D[:] = maketransmat(gpols[:,:,t],grid) * D[:];
    end
    

    return agg
end


"""numerical_jac compute a column of the Jacobians numerically
arguments:
ss steady state dictionary
T Int  length of transition path
j Int in [1, length(Xss)] specifying element of X to perturb
t Int   horizon at which to perturb, t = 1 means contemporaneous shock

return:
J  T x nagg array giving IRFs of each aggregate (where aggregates are computed by aggregator function)
"""
function numerical_jac(ss,T,j,t,par,grid)
    
    eps = 1e-4;

    Xss = ss[:X];
    Xpath0 = ones(T,1) .* Xss;

    Xpath1 = copy(Xpath0);
    Xpath1[t,j] += eps;
    

    agg0 =  solvesim(ss,Xpath0,par,grid);
    agg1 =  solvesim(ss,Xpath1,par,grid);

    return (agg1 .- agg0)./eps
end


function _fn_backwards_(ss,T,j,par,grid)
    
    eps = 1e-4;

    X0 = ones(T,1) .* ss[:X]; # this is for the ghost run
    X1 = copy(X0); # this is the perturbed version
    X1[T,j] += eps;
    
    dYs = zeros(T,ss[:Nagg]);
    dD1s = zeros(length(ss[:D]),T);
    
    # Solve backwards in time to get the savings policy rules
    Va_p_0 = copy(ss[:Va]); 
    Va_p_1 = copy(ss[:Va]);
    for t = T:-1:1
        Va_p_0,g0,c00 = EGMStepBack(Va_p_0,X0[t,:],par,grid);
        Va_p_1,g1,c1 = EGMStepBack(Va_p_1,X1[t,:],par,grid);
        
        # dYs[t,:] .= (aggregator(g1,c1,ss[:D]) .- aggregator(g0,c00,ss[:D])) / eps;
        dYs[t,:] .= (g1[:] .- g0[:])' * ss[:D] / eps;
        dD1s[:,t] = (maketransmat(g1,grid)  .- maketransmat(g0,grid) )* ss[:D] / eps; # ABRS equation 21: dD_1^s = dM_0^s D̄
    end
    return dYs, dD1s
end


function _fn_price_j_(ss,T,j,par,grid)

    dYs, dD1s = _fn_backwards_(ss,T,j,par,grid);

    # ABRS equation 23: E_t = M'^t y_ss
    M = maketransmat(ss[:g],grid);
    demean = x -> x .- sum(x)/length(x);

    ssaggvals = aggregator(ss[:g],ss[:c],1);  # we use 1 instead of D so we just get back the ss values associated with each grid point

    J = zeros(T,T,ss[:Nagg]);
    Es = zeros(length(ss[:D]),T-1);
    for i = 1:ss[:Nagg]

        Es[:].= 0.;
        if ss[:Nagg] == 1
            Es[:,1] .= demean(ssaggvals);
        else
            Es[:,1] .= demean(ssaggvals[i]);
        end
        for t = 2:T-1
            Es[:,t] = demean(M' * Es[:,t-1]);
        end



        # ABRS eq 24 (my timing is slightly different--my 1 is their 0)
        F = zeros(T,T);
        for t = 1:T
            for s = 1:T
                if t == 1
                    F[t,s] = dYs[T-(s-1),i]
                elseif t >= 2
                    F[t,s] = Es[:,t-1]' * dD1s[:,T-(s-1)]
                end
            end
        end
        J[:,:,i] = copy(F);
        for t = 2:T
            J[2:end, t,i] += J[1:end-1, t - 1,i]
        end
    end

    return J
end

function fakenews(ss,T,par,grid)
    Nx = length(ss[:X])
    J = zeros(T,T,ss[:Nagg],Nx);
    for j = 1:Nx
        J[:,:,:,j] = _fn_price_j_(ss,T,j,par,grid);
    end
    return J
end
