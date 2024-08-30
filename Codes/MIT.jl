using Plots
include("HA_fcns.jl")
include("HA_setup.jl")


compute_ss!(ss,par,grid);




"""Iterate over paths of capital to find an equilibrium"""
function MITequilibrium(ss,par,grid,Kguess,Zpath;maxit=100)
    
    

    for it =1:maxit
    
        savings = solvesim(ss,getWR(par,[ss[:K];Kguess[1:end-1]],Z=Zpath),par,grid);
    
       
        dif = maximum(abs.(Kguess .-savings))
        println("Update = $dif")
        if dif < 1e-6
            return Kguess
        else
            Kguess .+= 0.2*(savings .- Kguess)
        end
    end
    println("Did not converge")
end




T = 200
Zpath = exp.( 0.03 * 0.9.^(0:T-1)  );
Kguess = ss[:K]*ones(T)

eqmK = MITequilibrium(ss,par,grid,Kguess,Zpath)
Klag = [ss[:K];eqmK[1:end-1]];
println("K for t = 1 to 5: $(eqmK[1:5])")
X =getWR(par,Klag,Z=Zpath);
W = X[:,1];
R = X[:,2];
Lbar = 1.0;
Y = Zpath .* Klag.^par.α .* Lbar^(1-par.α);
    


function myplot(x,xss,lab)
    p = plot(0:60,x[1:61],title=lab,width=3,titlefontsize=12,xlabel="Quarter",label="");
    hline!([xss],label="",color=:black,width=2,linestyle=:dash)
    return p
end
series = [Zpath, eqmK, Y, W, R];
ssvals = [1,ss[:K],ss[:K]^par.α*Lbar^(1-par.α),getWR(par,ss[:K])...];
labs = ["Z", "K", "Y", "W","R"];
plts = [myplot(x...) for x in zip(series,ssvals,labs)];
l = @layout [a b; c d e]
plot(plts...,layout =l)
# savefig("Figures/MIT.png")