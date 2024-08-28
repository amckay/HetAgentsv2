using Plots
using Parameters


include("HA_fcns.jl")
include("HA_setup.jl")

# solve 
Va, g, c = SolveEGM(ss[:c],[W;R],par,grid);

# plot
display(plot(grid.a,c,title="Consumption Policy Rules",label=["Low endowment"  "High endowment"]))
ā = 100;
display(plot(grid.a[1:ā],[g[1:ā,:] agrid[1:ā]],title="Savings Policy Rules",label=["Low endowment"  "High endowment" "45-degree"],linestyles=[:solid :solid :dash],colors=[:red :blue  :black] ))
