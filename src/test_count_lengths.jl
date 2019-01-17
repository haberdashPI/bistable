using Pkg; Pkg.activate(joinpath(@__DIR__,".."))
include("count_lengths.jl")

datadir = joinpath(homedir(),"work","dlittle",
                   "bistable_individual_sensitive_W","run_2019-01-02")
# datadir = joinpath(@__DIR__,"..","data","count_lengths",
#                    "bistable_individual_levels","run_2018-11-15")

count_lengths(
  datadir=joinpath(datadir,"data"),
  logfile=joinpath(datadir,"logs","test.log"),
  first_index=1,last_index=1,
  sim_repeat=2,
  git_hash="UNKNOWN",
  params=joinpath(datadir,"params.jld2"),
  progressbar=false
)

