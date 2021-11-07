# Author: Raman Chhina
# Date: 2nd November 2021
# Purpose: Understanding Sorkin(2019)

using Distributions, LinearAlgebra, StatsBase, DataFrames, Plots, Chain, CategoricalArrays, FixedEffectModels, RegressionTables, LightGraphs
using Statistics, Random, SparseArrays, GLM

include("functions.jl")


mutable struct parameters
    α_sd::Float64   # Standard deviation of worker effects
    ψ_sd::Float64   # Standard deviation of firm effects
    V_sd:: Float64  # Standard deviation of Firm amenities
    ρ_V_Ψ::Float64  # Correlation between firm effect and firm amenity
    csort::Float64  # Sorting effect
    csig::Float64   # Cross-sectional standard deviation
    cnetw::Float64  # Network effects parameter
    w_sigma::Float64
    vsig::Float64
end

mutable struct types
    nk::Int     # approximate each distribution with some points of support
    nl::Int     # Let's assume moving probability is fixed
    λ::Float64  # Probability of switching
    γ::Float64  # Probability of death
    nt::Int     # numer of time periods
    ni::Int     # number of individuals
end


p = parameters(1.0,1.0,1.2,0.4,0.5,0.5,0.0,0.2,0.5)
q = types(30,10,0.1,0.05,10,1000)

α, ψ, V, G,H, B = transition_m(p::parameters, q::types)
df = simulate_panel(p, q, α, ψ, V, G, H, B)

df = lag_firms(df,10);
cc = connected_set(df)
df_connected = filter(y->y[1,:j] in cc[1],groupby(df, :i)) 
df_for_estimation = combine(df_connected, :i => same => :i, :j => same => :j, :t => same => :t, :lw => same => :lw, :α => same => :α, :ψ => same => :ψ)





