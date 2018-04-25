using Feather
include("count_lengths.jl")

parameter_index = length(ARGS) > 0 ? ARGS[1] : 1
dir = joinpath("..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

params = load("params.jld2")[:df]
N = 100
methods = Dict(
  :peaks => x -> source_count_by_peaks(x,window=1s,delta=0.25s,buildup=1s),
  :threshold => x -> source_count_by_threshold(x,window=1s,cuttoff=2cycoct,
                                               buildup=1s)
)
@assert log10(length(params)) < 8
count_lengths(@sprintf("params%08d",parameter_index),dir,N,methods,
              params[parameter_index])
