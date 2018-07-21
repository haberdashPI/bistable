using FileIO
using JLD2
using DataFrames
using Feather

include("count_lengths.jl")
params = Feather.read(params,transforms = Dict{String,Function}(
    "τ_x" => x -> totime.(x),
    "τ_σ" => x -> totime.(x),
    "τ_a" => x -> totime.(x),
    "τ_m" => x -> totime.(x),
    "delta_t" => x -> totime.(x),
    "standard_f" => x -> tofreq.(x),
    "condition" => x -> Symbol.(x)
  ))

# dir = joinpath("..","..","..","data","count_lengths")
dir = joinpath("work","dlittle","bistable_threshold_001","data")
all_rows = []
for_count_lengths(dir) do rows
  push!(all_rows,DataFrame(pindex = map(x -> x.pindex,rows),
                           length = map(x -> x.length,rows),
                           stimulus = map(x -> x.stimulus,rows),
                           created = map(x -> x.created,rows)))
end
df = vcat(all_rows...)

Feather.write(joinpath("..","..","..","data","count_lengths",
                       "scale_percept_lengths_$(Date(now())).feather"),df)

