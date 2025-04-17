using Plots
include("HA_fcns.jl")
include("HA_setup.jl")


compute_ss!(ss,par,grid);



# make the classic Aiyagari (1994) figure
Kgrid = KCompleteMarkets .+ collect(0:0.25:1.25) * (ss[:K] - KCompleteMarkets);
Ksupply = [K+checkK!(K,par,grid) for K in Kgrid];
Rgrid = [getWR(par,K)[2] for K in Kgrid];
plot(Ksupply,Rgrid,label="Capital supply",ylabel="R",xlabel="K",width=3);
p = plot!([Kgrid;48],[Rgrid;getWR(par,48)[2]],xlim=[40;60], label="Capital demand",legend=:bottomright,width=3)
display(p)
savefig("../Figures/AiyagariSD.png")

# plot the distribution of savings
normfactor = [1.0;diff(grid.a)];
density = reshape(ss[:D],length(grid.a),:) ./normfactor;
p = plot(grid.a,density,xlabel="Assets (before interest)",ylabel="Density", label=["Low endow" "High endow"],width=3,xlim=[-3;250])
display(p)
