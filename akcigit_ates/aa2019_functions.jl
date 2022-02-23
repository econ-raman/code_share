
# Function to calculate innovation rates given the Value Functions 

function calculate_innovation_rates(par::params, V::value_functions, ω::Float64)
    @unpack s, α, α̃, γ, γ̃, ϕ, ϕ̃, mbar = par
    @unpack v_l, v_f, v_nn = V

    x_l = zeros(mbar)
    x_f = zeros(mbar)
    xnn = 0
    x_e = zeros(mbar)
    x_enn = 0

    for m in 1:mbar-1
        x_l[m] = ((v_l[m+1] - v_l[m]) / ((1-s) * α * ω))^((1)/(γ-1))
    end

    for m in 2:mbar
        x_f[m] = ((max((ϕ*v_nn + (1- ϕ)*v_f[m-1] - v_f[m]),0) / ((1-s) * α * ω)))^((1)/(γ-1))
    end
    x_f[1] = (max((ϕ*v_nn + (1- ϕ)*v_nn - v_f[1]),0) / ((1-s) * α * ω))^((1)/(γ-1))

    xnn = ((v_l[1] - v_nn) / ((1-s) * α * ω))^((1)/(γ-1))

    for m in 2:mbar
        x_e[m] = ((ϕ̃ *v_nn + (1- ϕ̃ )*v_f[m-1] ) / ( α̃ * ω))^((1)/(γ̃ -1))
    end
    x_e[1] = ((ϕ̃ *v_nn + (1- ϕ̃ )*v_nn ) / ( α̃ * ω))^((1)/(γ̃ -1))


    x_enn = ((v_l[1] ) / ( α̃ * ω))^((1)/(γ̃ -1))

    return innovation_rates(x_l, x_f, xnn, x_e, x_enn)
end

# Functions to get the A's for value functions given the X's
function A_l(par::params, X::innovation_rates, V::value_functions, δ::Float64, ω::Float64)
    @unpack τ, λ, s, α, γ, ϕ, ϕ̃, mbar = par
    @unpack x_l, x_f, xnn, x_e, x_enn = X
    @unpack v_l, v_nn = V
    Al = zeros(mbar)
    for m in 2:mbar-1
        profit = (1-τ) * (1 - (1/ λ^m)) 
        continuation_1 = - ((1-s) * α * x_l[m]^(γ) * ω * (1/ γ)) +  x_l[m] * (v_l[m+1] - v_l[m])
        continuation_2 = ( ϕ * x_f[m] + δ + ϕ̃ * x_e[m] ) * (v_nn - v_l[m]) 
        continuation_3 =   ( (1-ϕ) * x_f[m] + (1-ϕ̃) * x_e[m] ) * (v_l[m-1] - v_l[m]) 
        Al[m] = profit + continuation_1 + continuation_2 + continuation_3 
    end

        profit = (1-τ) * (1 - (1/ λ)) 
        continuation_1 = - ((1-s) * α * x_l[1]^(γ) * ω * (1/ γ)) +  x_l[1] * (v_l[2] - v_l[1])
        continuation_2 = ( ϕ * x_f[1] + δ + ϕ̃ * x_e[1] ) * (v_nn - v_l[1]) 
        continuation_3 =   ( (1-ϕ) * x_f[1] + (1-ϕ̃) * x_e[1] ) * (v_nn - v_l[1]) 
        Al[1] = profit + continuation_1 + continuation_2 + continuation_3 

        profit = (1-τ) * (1 - (1/ λ)^mbar) 
        continuation_1 = - ((1-s) * α * x_l[mbar]^(γ) * ω * (1/ γ)) 
        continuation_2 = ( ϕ * x_f[mbar] + δ + ϕ̃ * x_e[mbar] ) * (v_nn - v_l[mbar]) 
        continuation_3 =   ( (1-ϕ) * x_f[mbar] + (1-ϕ̃) * x_e[mbar] ) * (v_l[mbar-1] - v_l[mbar]) 
        Al[mbar] = profit + continuation_1 + continuation_2 + continuation_3 
    return Al
end

function A_f(par::params, X::innovation_rates, V::value_functions, δ::Float64, ω::Float64)
    @unpack τ, λ, s, α, γ, ϕ, ϕ̃, mbar = par
    @unpack x_l, x_f, xnn, x_e, x_enn = X
    @unpack v_f, v_nn = V

    Af = zeros(mbar)
    for m in 2:mbar-1
        continuation_1 = - ((1-s) * α * x_f[m]^(γ) * ω * (1/ γ)) +  x_f[m] * (ϕ* v_nn + (1- ϕ) * v_f[m-1] - v_f[m])
        continuation_2 = δ * ( v_nn - v_f[m]) 
        continuation_3 =   x_l[m] * (v_f[m+1] - v_f[m])
        continuation_4 = x_e[m] *( 0 - v_f[m]) 
        Af[m] =   continuation_1 + continuation_2 + continuation_3 + continuation_4
    end

    continuation_1 = - ((1-s) * α * x_f[1]^(γ) * ω * (1/ γ)) +  x_f[1] * (ϕ* v_nn + (1- ϕ) * v_nn - v_f[1])
    continuation_2 = δ * ( v_nn - v_f[1]) 
    continuation_3 =   x_l[1] * (v_f[2] - v_f[1])
    continuation_4 = x_e[1] *( 0 - v_f[1]) 
    Af[1] =   continuation_1 + continuation_2 + continuation_3 + continuation_4

    continuation_1 = - ((1-s) * α * x_f[mbar]^(γ) * ω * (1/ γ)) +  x_f[mbar] * (ϕ* v_nn + (1- ϕ) * v_f[mbar-1] - v_f[mbar])
    continuation_2 = δ * ( v_nn - v_f[mbar]) 
    continuation_4 = x_e[mbar] *( 0 - v_f[mbar]) 
    Af[mbar] =   continuation_1 + continuation_2  + continuation_4

    return Af
end

function A_nn(par::params, X::innovation_rates, V::value_functions, ω::Float64)
    @unpack τ, λ, s, α, γ, ϕ, ϕ̃, mbar = par
    @unpack x_l, x_f, xnn, x_e, x_enn = X
    @unpack v_f, v_l, v_nn = V

    Ann = zeros(mbar)
    
        continuation_1 = - ((1-s) * α * xnn^(γ) * ω * (1/ γ)) +  xnn * (v_l[1]  - v_nn)
        continuation_2 = xnn * ( v_f[1] - v_nn) 
        continuation_3 =   x_enn* ( 0.5 * v_f[1] - v_nn)
        Ann = continuation_1 + continuation_2 + continuation_3 
    

    return Ann
end

# Functions to get the KFE A's
function KFE(par::params, X::innovation_rates, δ::Float64, μ::Array{Float64,1}, μ_nn::Float64)
    @unpack ϕ, ϕ̃ = par
    @unpack x_l, x_f, xnn, x_enn, x_e = X
   
    A_μ = zeros(mbar)
    A_μ[1] = μ_nn * (2 * xnn + x_enn) + μ[2] * ((1-ϕ)* x_f[2] + (1- ϕ̃)* x_e[2]) - μ[1] * (x_l[1] + x_f[1] + δ + x_e[1])
    for m = 2:mbar-1
        A_μ[m] = μ[m-1] * (x_l[m-1]) + μ[m+1] * ((1-ϕ)* x_f[m+1] + (1- ϕ̃)* x_e[m+1]) - μ[m] * (x_l[m] + x_f[m] + δ + x_e[m])   
    end
    A_μ[mbar] = μ[mbar-1] * x_l[mbar-1] - μ[mbar] *(x_f[mbar] + δ + x_e[mbar])
    A_μ_nn = μ[1] * ((1-ϕ) * x_f[1] + (1-ϕ̃)* x_e[1]) + sum( μ .* (ϕ .* x_f .+ δ .+ ϕ̃ .* x_e)) - μ_nn .* (2 * xnn + x_enn)

    return A_μ,A_μ_nn
end

# Function to solve the whole model given the value of δ
function solve_model(δ::Float64, par::params, ω_guess::Float64)

    ω = copy(ω_guess)
    dt = 1/50
    @unpack ρ = par 

    μ = (1/(mbar+1)) * ones(mbar) # Initial guess for the distribution 
    μ_nn = (1/(mbar))    # Initial guess for the neck to neck distribution
    V = value_functions(zeros(mbar), zeros(mbar),0) # Initial guess for the value functions 
    X = calculate_innovation_rates(par, V, ω)


    for p = 1:100 #loop for wage convergence

        V = value_functions(zeros(mbar), zeros(mbar),0) # Initial guess for the value functions 

        for t = 1:10000 #loop for VFI
            # First solve for the innovation rates
            X = calculate_innovation_rates(par, V, ω)

            # Given these innovation rates calculate updated value functions
            v_f_new = V.v_f .- (ρ .* V.v_f .- A_f(par,X,V, δ, ω)) .* dt
            v_l_new = V.v_l .- (ρ .* V.v_l .- A_l(par,X,V, δ, ω)) .* dt
            v_nn_new = V.v_nn - (ρ .* V.v_nn .- A_nn(par,X,V, ω)) .* dt

            # Check convergence 
            
            diff = maximum(abs.([v_f_new;v_l_new; v_nn_new] - [V.v_f;V.v_l;V.v_nn]))
            if diff < 10^-8
                break
            end

            #update value functions 
            V = value_functions(v_l_new, v_f_new, v_nn_new)
        end

        # Given this X solve for the stationary distribution 
            μ = (1/(mbar+1)) * ones(mbar)
            μ_nn = (1/(mbar))

            for i in 1:10000 #Loop for distribution convergence
                A_μ,A_μ_nn = KFE(par, X, δ, μ, μ_nn)
                μ_new = μ .+ A_μ .* dt
                #μ_nn_new = μ_nn + A_μ_nn .* dt
                diff = maximum(abs.(μ_new - μ))
                #diff = maximum(abs.([μ_new;μ_nn_new] - [μ;μ_nn]))
                if diff < 0.000001
                    break
                end
                μ = copy(μ_new)
                μ_nn = 1 - sum(μ)

            end

        @unpack λ, α, γ, α̃, γ̃ = par
        
        h_l = α .* X.x_l.^γ ./ γ
        h_f = α .* X.x_f.^γ ./ γ
        h_e = α .* X.x_e.^γ̃ ./ γ̃

        hl_0 =  α * ((X.xnn)^γ) / γ;
        hf_0 =  α * ((X.xnn)^γ) / γ
        he_0 = α̃ * ((X.x_enn)^γ̃) / γ̃

        temp0 = μ_nn * (hl_0 + hf_0 + he_0)
        den = 1 - sum(μ .* (h_l .+ h_f .+ h_e)) - temp0
        num = sum([μ[m] / λ^m for m in 1:mbar]) + μ_nn
        ω_that_clears_market = num /den


        diff = abs(ω_that_clears_market - ω)
        if diff < 1e-4
            break
        end
            ω = 0.25 * ω_that_clears_market + 0.75 * ω
        println(diff)
    end
    @unpack λ = par
    g = log(λ) * (μ_nn *(2* X.xnn + X.x_enn ) + sum(μ .* X.x_l))
    return g, V, X
end