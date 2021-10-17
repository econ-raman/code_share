function basic_checks(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)

    @assert size(a,1) > 2 "a should have atleast two observations"
    @assert size(a,2) == 2 "a should have atleast two columns"
    @assert all(a .> 0) "all the entries of a should be strictly positive"

    A = a[:,1]./a[:,2];

    for i in 1:(size(A,1) - 1)
        @assert (A[i] - A[i+1] > 0) "A is not monotonically decreasing"
    end

    @assert (size(b,1) == size(a,1)) "b should have the same length as a"
    @assert all(b .> 0) "all the elements of b should be stictly positive"
    @assert round(sum(b),digits = 2) == 1.0 "elements of b should sum to 1"

    @assert g <= 1 "g should be a scalar less than equal to 1"
    @assert g > 0 "g should be strictly positive"

end

    

function DFS1977fig1(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1}, g::Float64)
    basic_checks(a,b,L,g)
    A = a[:,1]./a[:,2];
    N = length(A);
    fig1 = plot()
    z = collect(range(0, stop = 1, length = N));
    fig1 = plot!(z,A, xlabel = "z", label = "A(z)", linewidth = 3, palette = :grays)

    ν(z̃) = sum(b[1:Int(floor(z̃*N))])
    B(z̃) = (ν(z̃)/(1-ν(z̃))) * (L[1] / L[2]);
    z_capped = z[1:Int(floor(0.9*N))]
    fig1 = plot!(z_capped, B.(z_capped), label = "B(z;L*/L)", linewidth = 3, palette = :grays, legend = :topleft)
    return fig1
end

function DFS1977solver_g1(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)

  basic_checks(a,b,L,g)
  A = a[:,1]./a[:,2];
  N = length(A);
  z = collect(range(0, stop = 1, length = N));

  ν(z̃) = sum(b[1:Int(floor(z̃*N))])
  B(z̃) = (ν(z̃)/(1-ν(z̃))) * (L[1] / L[2]);

  # simple grid search Solution
  
    z̄ = 0
    F(i) = B(z[i]) - A[i]
  for i in 2:N
        if ((F(i) * F(i-1)) < 0)
            z̄  = i
            break
        end
  end
  ω̄ = (A[z̄] + A[z̄-1]) / 2
  return ω̄, z[z̄]

end

function z_bar(g::Float64,ω::Float64, z::Array{Float64,1}, A::Array{Float64,1}) 
    Az = LinearInterpolation(z, A)
    f(z̄) =  Az(z̄) - (g*ω);
    z̄ = find_zero(f, (0,1))
    return z̄
end

function z_bar_star(g::Float64,ω::Float64, z::Array{Float64,1}, A::Array{Float64,1}) 
    Az = LinearInterpolation(z, A)
    f(z̄) =  Az(z̄) - (ω/g);
    z̄ = find_zero(f, (0,1))
    return z̄
end

function DFS1977solver_any_g(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)

    basic_checks(a,b,L,g)
    A = a[:,1]./a[:,2];
    N = length(A);
    z = collect(range(0, stop = 1, length = N));
    
    A_1 = A ./ g;

         λ(g, ω) = sum(b[1:Int(floor(z_bar(g,ω,z,A)  *N))])
    λ_star(g, ω) = sum(b[Int(floor(z_bar(g,ω,z,A)  *N)):end])

  
    # simple grid search Solution

     ω̄_grid = collect(range(A_1[1], stop = A_1[end], length = 1000))
     F(ω̄) =  ω̄ - (( (1 - λ_star(g, ω̄)) / (1 - λ(g, ω̄)) ) * (L[1] / L[2]))
     ω_bar = 0.0
    for i in 1:99
        if (F(ω̄_grid[end+1-i]) * F(ω̄_grid[end-i]) < 0)
            break
        end
        ω_bar = (ω̄_grid[end+1-i] +  ω̄_grid[end-i]) / 2;
    end
    

  
    z̄ = z_bar(g, ω_bar, z, A)
    z̄_star = z_bar_star(g, ω_bar, z, A) 

    return ω_bar, z̄, z̄_star
  
  end

function DFS1977solver(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)
    if g == 1.0
        DFS1977solver_g1(a,b,L,g)
    else
        DFS1977solver_any_g(a,b,L,g)
    end
end
