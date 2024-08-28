using Plots
include("HA_fcns.jl")
include("HA_setup.jl")


compute_ss!(ss);






T = 200;
print("Starting fake news...")
J = fakenews(ss,T)
println("done")


print("Starting numerical Jacobian columns...")
nm_W = numerical_jac(ss,T,1,10);
nm_R = numerical_jac(ss,T,2,10);
println("done")

p_W = plot([J[1:100,10,1,1] nm_W[1:100]], labels=["Fake news" "Numerical"], width=3, title="W", titlefontsize=12,linestyles=[:solid :dash])
p_R = plot([J[1:100,10,1,2] nm_R[1:100]], labels=["Fake news" "Numerical"], width=3, title="R", titlefontsize=12,linestyles=[:solid :dash])
plot(p_W, p_R)
savefig("Figures/SSJ.png")