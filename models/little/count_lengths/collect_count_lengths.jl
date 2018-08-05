using FileIO
using JLD2
using DataFrames
using Feather

include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
# dir = joinpath("work","dlittle","bistable_threshold_001","data")
all_rows = []
for_count_lengths(dir) do count_length
  push!(all_rows,DataFrame(pindex = count_length.pindex,
                           created = count_length.created,
                           ratio = count_length.ratio,
                           kind = "component"))
  push!(all_rows,DataFrame(pindex = count_length.pindex,
                           created = count_length.created,
                           ratio = count_length.bratio,
                           kind = "bandwidth"))
end
df = vcat(all_rows...)

Feather.write(joinpath("..","..","..","data","count_lengths",
                       "scale_percept_lengths_$(Date(now())).feather"),df)

