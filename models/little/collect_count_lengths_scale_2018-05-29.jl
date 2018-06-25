using FileIO
using JLD2
using DataFrames
using Feather

include("count_lengths.jl")
params = load("params.jld2","df")

dir = joinpath("..","..","data","count_lengths","data")
all_rows = []
for_count_lengths(dir) do rows
  if !any(x -> x.method ∉ 1:3,rows)
    push!(all_rows,DataFrame(pindex = map(x -> x.pindex,rows),
                             length = map(x -> x.length,rows),
                             stimulus = map(x -> x.stimulus,rows),
                             method = map(x -> method_labels[x.method],
                                          rows)))
  end
end
df = vcat(all_rows...)

Feather.write(joinpath("..","..","data","count_lengths",
                       "scale_percept_lengths_2018-05-29.feather"),df)

