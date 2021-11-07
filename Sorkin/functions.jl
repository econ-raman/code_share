# Author: Raman Chhina
# Date: 2nd November 2021
# Purpose: Functions for Understanding Sorkin(2019)

function transition_m(p::parameters, q::types)
    α_sd = p.α_sd; ψ_sd = p.ψ_sd; csort = p.csort; cnetw = p.cnetw; csig = p.csig; ρ_V_Ψ = p.ρ_V_Ψ; V_sd = p.V_sd; vsig = p.vsig;
    nk = q.nk; nl = q.nl; λ = q.λ; 

    ψ = quantile.(Normal(), (1:nk) / (nk + 1)) * ψ_sd
    α = quantile.(Normal(), (1:nl) / (nl + 1)) * α_sd
    V = quantile.(Normal(), (1:nk) / (nk + 1)) * V_sd

    # Given a value of Ψ we define a distribution for V_sd
    B = zeros(nk,nk)
    for k in 1:nk
        B[k,:] = pdf(Normal(0,vsig), V .- ρ_V_Ψ * ψ[k]) # Firm effect and amenity are correlated
        B[k,:] = B[k,:] ./ sum(B[k,:])
    end

    # Let's create type-specific transition matrices
    # We are going to use joint normals centered on different values
    G = zeros(nl, nk, nk)
    for l in 1:nl, k in 1:nk
        G[l, k, :] = pdf(Normal(0, csig), ψ .- cnetw * ψ[k] .- csort * α[l])
        G[l, k, :] = G[l, k, :] ./ sum(G[l, k, :])
    end

    # We then solve for the stationary distribution over psis for each alpha value
    # We apply a crude fixed point approach
    H = ones(nl, nk) ./ nk
    for l in 1:nl
        M = transpose(G[l, :, :])
        for i in 1:100
            H[l, :] = M * H[l, :]
        end
    end

  
    return α, ψ, V, G,H, B
end

function simulate_panel(p::parameters, q::types, α::Array{Float64,1}, ψ::Array{Float64,1}, V::Array{Float64,1}, G::Array{Float64, 3}, H::Array{Float64,2}, B::Array{Float64,2})
    α_sd = p.α_sd; ψ_sd = p.ψ_sd; csort = p.csort; cnetw = p.cnetw; csig = p.csig; γ = q.γ;
    nk = q.nk; nl = q.nl; λ = q.λ; ni = q.ni; nt = q.nt;
    
    # We simulate a balanced panel
    ll = zeros(Int64, ni, nt) # Worker type
    kk = zeros(Int64, ni, nt) # Firm type
    vv = zeros(Int64, ni, nt) # Amenity type
    spellcount = zeros(Int64, ni, nt) # Employment spell

    ii = zeros(Int64, ni, nt) # Worker Ids
    p = ni + 1; # counter for new workers who get added after a worker dies

    for i in 1:ni
        ii[i,:] .= i
        # We draw the worker type
        l = rand(1:nl)
        ll[i,:] .= l
        
        # At time 1, we draw from H
        kk[i,1] = sample(1:nk, Weights(H[l, :]))
        # Given firm type we draw amenity
        vv[i,1] = sample(1:nk, Weights(B[kk[i,1],:]))
        draw1 = rand()
        draw2 = rand()
        for t in 2:nt
            if ((draw1 < λ) & (draw2 > γ)) # Doesn't dies but gets a switch offer
                # draw the new firm and amenity
                k_new= sample(1:nk, Weights(G[ll[i,t], kk[i,t-1], :]))
                v_new = sample(1:nk, Weights(B[k_new,:]))
                # move with probability
                if rand() < ((exp(v_new))/(exp(v_new) + exp(vv[i,t-1])))
                    kk[i,t] = k_new
                    vv[i,t] = v_new
                    spellcount[i,t] = spellcount[i,t-1] + 1
                else
                    kk[i,t] = kk[i,t-1]
                    vv[i,t] = vv[i,t-1]
                    spellcount[i,t] = spellcount[i,t-1]
                end
            elseif ((draw1 > λ) & (draw2 > γ))
                kk[i,t] = kk[i,t-1]
                vv[i,t] = vv[i,t-1]
                spellcount[i,t] = spellcount[i,t-1]
            else
                #drop the worker and redraw the type, firm and amenity
                ii[i, t:end] .= p; p = p + 1; # give them a new id and increase the new worker counter
                l = rand(1:nl)
                ll[i,t:end] .= l
                kk[i,t] = sample(1:nk, Weights(H[ll[i,t], :]))
                vv[i,t] = sample(1:nk, Weights(B[kk[i,t],:]))
            end
        end
        
    end

    # Attach Firm ID to types
    # The final step is to assign identities to the firms. We are going to do this is a relatively simple way, by simply randomly assigning firm ids to spells.
    firms_per_type = 15
    jj = zeros(Int64, ni, nt) # Firm identifiers

    draw_firm_from_type(k) = sample(1:firms_per_type) + (k - 1) * firms_per_type

    for i in 1:ni
        
        # extract firm type
        k = kk[i,1]
        
        # We draw the firm (one of firms_per_type in given group)
        jj[i,1] = draw_firm_from_type(k)
        
        for t in 2:nt
            if spellcount[i,t] == spellcount[i,t-1]
                # We keep the firm the same
                jj[i,t] = jj[i,t-1]
            else
                # We draw a new firm
                k = kk[i,t]
                
                new_j = draw_firm_from_type(k)            
                # Make sure the new firm is actually new
                while new_j == jj[i,t-1]
                    new_j = draw_firm_from_type(k)
                end
                
                jj[i,t] = new_j
            end
        end
    end

    # Make sure firm ids are contiguous
    contiguous_ids = Dict( unique(jj) .=> 1:length(unique(jj))  )
    jj .= getindex.(Ref(contiguous_ids),jj);

    # Construct DataFrame
    #ii = repeat(1:ni,1,nt) - now I'm using the new definition with new workers
    tt = repeat((1:nt)',ni,1)
    df = DataFrame(i=ii[:], j=jj[:], l=ll[:], k=kk[:], α=α[ll[:]], ψ=ψ[kk[:]], V=V[vv[:]], t=tt[:], spell=spellcount[:]);
    return df
end

function lag_firms(df::DataFrame, T::Int64)
    lag_firm = zeros(Int64, q.ni*q.nt)
    sort!(df, [:i,:t])
    for i in 1:q.ni
        sub_df = df[df.i .== i,:]
        k = ((i-1)*T) + 1
        lag_firm[k] = sub_df.j[1]
        lag_firm[(k+1):(k+T-1)] = sub_df.j[1:end-1]
    end

    df[!,:lag_j] = lag_firm ;
    return df
end


function connected_set(df::DataFrame)
    nj = length(unique(df.j))
    adjacency_matrix = zeros(Int, nj, nj);
    
    for a in unique(df.j), b in unique(df.j)
      if (size((df[(df.j .== a) .& (df.lag_j .== b),:]),1) !== 0)
        adjacency_matrix[a,b] = 1; adjacency_matrix[b,a] = 1;
      end
  end

  g = SimpleGraph(adjacency_matrix);
  cc = connected_components(g)

return cc
end