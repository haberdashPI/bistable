using Feather
include("count_lengths.jl")
# NEXT STEP: find all using statements and install those packages
# on the MARCC cluster
parameter_index = length(ARGS) > 0 ? ARGS[1] : 1
dir = length(ARGS) > 1 ? ARGS[2] : joinpath("..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

params = load("params.jld2")["df"]
N = 25
M = 10
#  N = 1
#  M = 1

methods = Dict(
  # :peaks => x -> source_count_by_peaks(x,window=1s,delta=0.25s,buildup=1s),
  :threshold => x -> source_count_by_threshold(x,window=1s,delta=0.25s,
                                               cutoff=2cycoct,buildup=1s)
)

@assert log10(length(params)) < 8
count_lengths(@sprintf("params%08d",parameter_index),dir,N,M,methods,
              Dict(k => params[parameter_index,k] for k in names(params)))
