# Author: Raman Chhina
# Purpose: Replication of Akcigit and Ates (2019)
# Date Last Modified: 23rd February 2022

include("./aa2019_functions.jl") # All the functions are stored here
using Parameters, Roots

# Define a structure to hold the model parameters
@with_kw struct params
    ρ::Float64 = 0.05
    γ::Float64 = 1 / 0.35
    γ̃::Float64 = 1 / 0.35
    τ::Float64 = 0.3
    s::Float64 = 0.05
    α::Float64 = 7.179
    α̃::Float64 = 0.075
    λ::Float64 = 1.044
    ϕ::Float64 = 0.0423
    ϕ̃::Float64 = 0.0423
    mbar::Int64 = 100
end

# Structures to hold the valye functions and innovation rates

mutable struct value_functions 
    v_l::Array{Float64,1}  # Value function of the leader
    v_f::Array{Float64,1}  # Value function of the follower
    v_nn::Float64          # Value of the neck to neck firm 
end

mutable struct innovation_rates
    x_l::Array{Float64,1}  # Leader innovation rate
    x_f::Array{Float64,1}  # Follower innovation rate
    xnn::Float64           # Neck to neck incumbent innovation rate
    x_e::Array{Float64,1}  # Entrant Innovation rate
    x_enn::Float64         # Innovation rate of an entrant in the neck to neck sector
end



par = params() # Initialise the parameters 

@unpack mbar = par # Get mbar
μ = zeros(mbar)   # Mass of unleveled sectors
μ_nn = 0          # Neck to Neck sectors

#Define the grid on δ
δ_grid = collect(range(0.015, 5, length = 15))
δ = δ_grid[11]

# Initial guess for ω
ω_guess = 0.9

gg = zeros(15) # Holder for growth rates

# Solve the model for each δ and store the g
for j in 1:15
    gg[j], V, X = solve_model(δ_grid[j], par,0.9)
end

# Plot the final output
using Plots
plot(δ_grid,  gg, linewidth = 3, grid = :none, label = "", xlabel = "Delta", ylabel = "g")
savefig("plot.png")