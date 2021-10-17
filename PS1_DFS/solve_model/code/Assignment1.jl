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




fig1 = DFS1977fig1(a,b,L)
ω̄, z̄ = DFS1977solver(a,b,L,g)
fig1 = vline!([z̄], label = "z̄")
fig1 = hline!([ω̄], label = "ω̄")
savefig(fig1, "fig1.png")


a = readdlm("/home/rschhina/second_year/second_year/International_Trade/Assignments/PS1_DFS/solve_model/input/DFS1977_example_a.txt",  '\t', Float64)
b = vec(readdlm("/home/rschhina/second_year/second_year/International_Trade/Assignments/PS1_DFS/solve_model/input/DFS1977_example_b.txt",  '\t', Float64))


omega_bar, zbar = DFS1977solver(a,b,L,1.0)
println("Solution: z_bar_star = ", zbar*N, " z_bar = ", zbar*N, " omega_bar = ", omega_bar)
omega_bar, zbar,zbar_star = DFS1977solver_any_g(a,b,L,0.9)
println("Solution: z_bar_star = ", zbar_star*N, " z_bar = ", zbar*N, " omega_bar = ", omega_bar)


