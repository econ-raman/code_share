# Find all movers defined as those with max spell > 0
function find_movers(df::DataFrame)
    move_df = @chain df begin
       groupby(:i)
       transform(:spell => maximum) 
       subset(:spell_maximum => x -> x .> 0)
    end

   return move_df 
end


# Get the next firm the mover is moving to
# by leading the firm column and extracting 
# just the firm and next firm as one observation
function find_firm_links(mover_df::DataFrame)
    firm_link_df = @chain mover_df begin
        sort([:i, :t])
        groupby(:i)
        transform(:j => lead => :j_next)
        transform([:j, :j_next] => .==)
        subset(:j_j_next_BroadcastFunction => x -> x .== false, skipmissing = true)
        select(:j, :j_next)
        unique()
    end
    return firm_link_df
end

# Iterate through our firm link df and create a matrix 
# of links. In R I would model.matrix(a ~ b) but idk how to 
# do that here
function create_adjacency_matrix(firm_link_df::DataFrame, df::DataFrame)
    adjacency_matrix = zeros(Int, maximum(df.j), maximum(df.j))
    for firm_a in unique(df.j), firm_b in unique(df.j)
        subset_a_df = firm_link_df[(firm_link_df.j .== firm_a) .& (firm_link_df.j_next .== firm_b), :]
        if size(subset_a_df)[1] != 0
            adjacency_matrix[firm_a, firm_b] = 1
            adjacency_matrix[firm_b, firm_a] = 1
        end
    end
    return adjacency_matrix
end

# %%
mover_df = find_movers(df)
firm_link_df = find_firm_links(mover_df)
adjacency_matrix = create_adjacency_matrix(firm_link_df, mover_df)
