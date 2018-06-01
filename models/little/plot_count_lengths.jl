using FileIO
using JLD2
using DataFrames
using Feather

include("count_lengths.jl")
params = load("params.jld2","df")

# start by just looking the mean and max length
# of the two percepts across the various
# parameters

dir = joinpath("..","..","data","count_lengths")
all_rows = []
for_count_lengths(dir) do rows
  if !any(row -> row.error_code > 0,rows)
    push!(all_rows,DataFrame(pindex = map(x -> x.pindex,rows),
                             length = map(x -> x.length,rows),
                             stimulus = map(x -> x.stimulus,rows),
                             created = map(x -> x.created,rows),
                             method = map(x -> method_labels[x.method],rows)))
  end
end
df = vcat(all_rows...)

Feather.write(joinpath("..","..","data","count_lengths","all_rows.feather"),df)

# r_params = DataFrame(
#   c_sigma = params[:c_σ],
#   tau_sigma = params[:τ_σ],
#   c_m = params[:c_m],
#   tau_m = params[:τ_m],
#   c_a = params[:c_a],
#   tau_a = params[:τ_a],
#   W_m_sig = params[:W_m_σ],
#   tau_x = params[:τ_x],
#   c_x = params[:c_x],
#   c_n = params[:c_n],
#   τ_n = params[:τ_n],
# )
# Feather.write(joinpath("..","..","data","count_lengths","params.feather"),r_params)


