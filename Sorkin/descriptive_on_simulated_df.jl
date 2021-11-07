# Author: Raman Chhina
# Date: 2nd November 2021
# Purpose: Understanding Sorkin(2019): Decriptive analysis on the simulated data

using Distributions, LinearAlgebra, StatsBase, DataFrames, Plots, Chain, CategoricalArrays, FixedEffectModels, RegressionTables, LightGraphs
using Statistics, Random, SparseArrays, GLM
default( palette = :grays, label = "", size = (800,300))


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

include("functions.jl")

p = parameters(1.0,1.0,1.2,0.4,0.5,0.5,0.0,0.2,0.5)
q = types(30,10,0.1,0.05,10,1000)

α, ψ, V, G,H, B = transition_m(p::parameters, q::types)
df = simulate_panel(p, q, α, ψ, V, G, H, B)


# Calculate correlation between Ψ and V
cor(df.ψ, df.V) # I get a value of around 0.4 as expected by the parameterisation. Good!

# Number of worker deaths
size(unique(df.i),1) - q.ni # I get a value of around 500 as there are 1000 individuals and 10 time periods this is Good

# correlation between worker and firm effects
cor(df.ψ, df.α)

# ------------------------------ Number of per movers per firm -----------------------------------#
p.vsig = 1
p.V_sd = 1
q.λ = 0.3
df = simulate_panel(p, q, α, ψ, V, G, H, B)

df_movers = combine(groupby(df,:i), :j => length ∘ unique)
rename!(df_movers, :j_ComposedFunction => :lifetime_firms)
histogram(df_movers.lifetime_firms)
savefig("figures/lifetime_firms.png")


df = join(df, df_movers, on = :i, kind = :left);

movers_per_firm = filter(row -> row.lifetime_firms > 1, df)
movers_per_firm  = combine(groupby(df,[:j, :t]), :i => length ∘ unique)
rename!(movers_per_firm, :i_ComposedFunction => :movers_in_this_firm)


cross_section = movers_per_firm[movers_per_firm.t .== 10,:]
histogram(cross_section.movers_in_this_firm)
mean_movers = round(mean(cross_section.movers_in_this_firm), digits = 3)
vline!([mean_movers], linewidth = 3, label = "mean = $mean_movers", palette = :grays, legend = :topleft)
title!("Distribution of movers per firm in the cross section at t = 5")
savefig("figures/movers_per_firm.png")


