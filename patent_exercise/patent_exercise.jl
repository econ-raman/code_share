using DataFrames, CSV, Dates, Statistics, ShiftedArrays, Plots
ENV["GKSwstype"] = "100";
Plots.default(legend = :outerright, palette = :grays, linewidth = 3, linecolor = :black, ls = :auto, grid = :none, size = (800,400));

# Read Data
assignee = CSV.read("/home/rschhina/income_inequality/patent/assignee.tsv", delim = "    ")
inventor = CSV.read("/home/rschhina/income_inequality/patent/inventor.tsv", header = ["id", "name_first", "name_last", "male_flag", "attribution_status"], datarow = 2)
patent_assignee = CSV.read("/home/rschhina/income_inequality/patent/patent_assignee.tsv", header=["patent_id", "assignee_id", "location_id"], datarow=2)
patent_inventor = CSV.read("/home/rschhina/income_inequality/patent/patent_inventor.tsv", header=["patent_id", "inventor_id", "location_id"], datarow=2)
patent = CSV.read("/home/rschhina/income_inequality/patent/patent.tsv", header = ["id", "type", "number", "country", "date", "abstract", "title", "kind", "num_claims", "filename", "withdrawn"], datarow = 2)

# Merge the date info of the patent with `patent-assignee` crosswalk
patent_date = select(patent, [:id, :number, :date])
rename!(patent_date, :id => :patent_id)
patent_date.year = Dates.year.(patent_date.date)
patent_assignee_date = join(patent_date, patent_assignee, on = :patent_id, kind = :left)

# Remove Duplicates
patent_assignee_date_unique = unique(patent_assignee_date)

# Count the total number of patents for each firm between 1980 and 2018
patent_firms = combine(groupby(patent_assignee_date_unique, [:assignee_id]), nrow)
rename!(patent_firms, :x1 => :total_patents)
sort!(patent_firms, :total_patents, rev = :true)

# Define the top 1 percent firms 
q99 = quantile(patent_firms.total_patents, 0.99)
patent_firms.top_firm = ifelse.(patent_firms.total_patents .> q99, 1, 0) 

# Each year how many inventors move from the top firms to laggard firms 

# Find movers

# Steps
# For each patent get the assignee and the inventor together 
# Find inventors who are listed with more than one assignee

rename!(inventor, :id => :inventor_id)

patent_inventor_firms = join(patent_assignee_date_unique, patent_inventor, on = :patent_id, kind = :left, makeunique = true)
select!(patent_inventor_firms, [:patent_id, :date, :year, :assignee_id, :inventor_id])

same(x) = x
inventor_firms = combine(groupby(patent_inventor_firms, [:inventor_id, :date]), :assignee_id => same )
unique!(inventor_firms)
rename!(inventor_firms, :assignee_id_same => :assignee_id)
inventor_firms = join(inventor_firms, patent_firms, on = :assignee_id, kind = :left)
sort!(inventor_firms, [:inventor_id, :date])

# Remove Missing

inventor_firms_clean = inventor_firms[inventor_firms.inventor_id .!== missing,:]
inventor_firms_clean_1 = inventor_firms_clean[inventor_firms_clean.assignee_id .!== missing,:]

# Find Movers
inventor_firms_clean_1.year = Dates.year.(inventor_firms_clean_1.date)


# Find the lag firms for each inventor
sort!(inventor_firms_clean_1, [:inventor_id, :date])
l1(x) = lag(x, 1)
inventor_firms_clean_2 = by(inventor_firms_clean_1, :inventor_id,  :assignee_id => l1 ,
:top_firm => l1  , 
 :total_patents => l1  , 
 :assignee_id => same,
  :total_patents => same,
  :date => same, 
  :top_firm => same )

# Keep only the observations where lag assignee id is not missing. This automatically removes the inventors 
# who are with only one assignee. And in other cases removes the first time the inventor is seen in the data. 
inventor_firms_clean_atleat_two = inventor_firms_clean_2[inventor_firms_clean_2.assignee_id_l1 .!== missing,:]



function moves(j, lag_j)
    if (j == lag_j)
        move = "No"
    else
        move = "Yes"
    end
    return move
end

rename!(inventor_firms_clean_atleat_two,
:assignee_id_same => :assignee_id, 
:total_patents_same => :total_patents, 
:date_same => :date, 
:top_firm_same => :top_firm)

# Find the time periods in which the inventor moved
inventor_firms_clean_atleat_two[!, :move] = moves.(inventor_firms_clean_atleat_two.assignee_id, inventor_firms_clean_atleat_two.assignee_id_l1)
mover_times = select( inventor_firms_clean_atleat_two[inventor_firms_clean_atleat_two.move .== "Yes",:], [:inventor_id, :date, :assignee_id, :assignee_id_l1, :top_firm, :top_firm_l1, :total_patents, :total_patents_l1])
rename!(mover_times, :date => :move_date) ;

# Label move as 1 if it is from a top firm to laggard firm and 0 otherwise
function top_to_laggard(j, lag_j)
    if (j == 0) & (lag_j == 1)
        move = 1
    else
        move = 0
    end
    return move
end

mover_times.top_to_laggard = top_to_laggard.(mover_times.top_firm, mover_times.top_firm_l1)

mover_times.move_year = Dates.year.(mover_times.move_date)

# Calculate total moves in a year and fraction which is top to laggard.
total_laggard_moves = combine(groupby(mover_times, :move_year), :top_to_laggard => sum)
total_moves = combine(groupby(mover_times, :move_year), nrow)

sort!(total_laggard_moves, :move_year)
sort!(total_moves, :move_year)


plot(total_laggard_moves.move_year, total_laggard_moves.top_to_laggard_sum, label = "top to laggard" )
plot!(total_laggard_moves.move_year,  total_moves.x1, label = "total")
savefig("total_moves")

plot(total_laggard_moves.move_year, total_laggard_moves.top_to_laggard_sum ./ total_moves.x1, label = "top to laggard fraction")
savefig("fraction.png")

