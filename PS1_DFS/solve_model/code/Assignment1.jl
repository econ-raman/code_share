##Dingel, DFS 1977, example of function call
##Set up Julia packages
import Pkg
Pkg.activate("../input/Project.toml")
Pkg.instantiate()
ENV["GKSwstype"] = "100";

using CSV, DelimitedFiles, LaTeXStrings, Interpolations, Random, Plots, Roots  #Load packages

include("DFS1977functions.jl")


Random.seed!(1234) # set seed
N = 150; # no. of varieties
a =hcat(ones(N,1),0.25*ones(N,1) + rand(N,1));
a[:,2] = sort(a[:,2])
b = 0.5*rand(N) + 0.5*ones(N)
b = b / sum(b)
L = ones(2)
g = 1.0;


a = readdlm("../input/DFS1977_example_a.txt",  '\t', Float64)
b = vec(readdlm("../input/DFS1977_example_b.txt",  '\t', Float64))

# Comparison to results in the Assignment PDF
omega_bar, zbar = DFS1977solver(a,b,L,1.0)
println("Solution: z_bar_star = ", zbar*N, " z_bar = ", zbar*N, " omega_bar = ", omega_bar)

# No trade cost equilibrium
fig1 = plot()
fig1 = DFS1977fig1(fig1,a,b,L,g)
ω̄, z̄ = DFS1977solver(a,b,L,g)
fig1 = vline!([z̄], label = "z̄", ls = :dash, linecolor = :black)
fig1 = hline!([ω̄], label = "ω̄", ls = :dot, linecolor = :black)
savefig(fig1, "../output/fig1.png")


# Comparative Statistics

# Increase in foreign size
L_new = [1.5;1]
fig2 = DFS1977fig1(fig1,a,b,L_new ,g)
ω̄, z̄ = DFS1977solver(a,b,L,g)
fig2 = vline!([z̄], label = "z̄'", ls = :dash, linecolor = :black)
fig2 = hline!([ω̄], label = "ω̄'", ls = :dot, linecolor = :black)
savefig(fig2, "../output/fig2.png")

# Uniform foreign technical progress
a_new = copy(a);
a_new[:,1] = 1.2 .* a[:,1]
fig2_1 = DFS1977fig1(fig1,a_new,b,L ,g)
ω̄, z̄ = DFS1977solver(a_new,b,L,g)
fig2_1 = vline!([z̄], label = "z̄'", ls = :dash, linecolor = :black)
fig2_1 = hline!([ω̄], label = "ω̄'", ls = :dot, linecolor = :black)
savefig(fig2_1, "../output/fig2_1.png")

# Introduction of trade costs
omega_bar, zbar,zbar_star = DFS1977solver_any_g(a,b,L,0.9)
println("Solution: z_bar_star = ", zbar_star*N, " z_bar = ", zbar*N, " omega_bar = ", omega_bar)

# Fig 3 of DFS 1977
fig3 = DFS_fig3(a,b,L,0.7)
savefig(fig3, "../output/fig3.png")

# Gains from trade when foreign experiences a uniform technical progress
home_gainf(y) = home_welfare(y,b,L,g) - home_welfare(a,b,L,g)
foreign_gainf(y) = foreign_welfare(y,b,L,g) - foreign_welfare(a,b,L,g)

a_increase = collect(range(1, stop = 3, step = 0.1))

home_gains = zeros(size(a_increase,1))
foreign_gains = zeros(size(a_increase,1))


for (i,j) in zip(a_increase, 1:length(a_increase))
    home_gains[j] = home_gainf([a[:,1].*i  a[:,2]])
    foreign_gains[j] = foreign_gainf([a[:,1].*i  a[:,2]])
end

w_plot = plot(a_increase, [ home_gains,  foreign_gains], linewidth = 3, ls = [:dash :dot], linecolor = :black,
label = ["Home gain" "Foreign gain"], legend = :topleft, xlabel = "Proportion increase in a*(z)",
ylabel = "Gain")
savefig(w_plot, "../output/welfare.png")

# Gains and trade_volume

current_vot = trade_volume(a,b,L,g)
println("The volume of trade is ", current_vot)

# Find another equilibrium (with same L and g) which has same vot and different gains 

home_gains = zeros(size(a_increase,1))
foreign_gains = zeros(size(a_increase,1))
vot = zeros(size(a_increase,1))


for (i,j) in zip(a_increase, 1:length(a_increase))
    home_gains[j] = home_gainf([a[:,1].*i  a[:,2]])
    foreign_gains[j] = foreign_gainf([a[:,1].*i  a[:,2]])
    vot[j] = trade_volume([a[:,1].*i  a[:,2]], b, L,g)
end

w_plot = plot(a_increase, [ home_gains,  foreign_gains, vot], linewidth = 3, ls = [:dash :dot :solid], linecolor = :black,
label = ["Home gain" "Foreign gain" "Volume of Trade"], legend = :topleft, xlabel = "Proportion increase in a*(z)",
ylabel = "Gain")
savefig(w_plot, "../output/welfare_vot.png")


