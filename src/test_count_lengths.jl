using Pkg; Pkg.activate(joinpath(@__DIR__,".."))
include("count_lengths.jl")

datadir = normpath(joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-06"))
count_lengths(
  datadir=joinpath(datadir,"data"),
  logfile=joinpath(datadir,"logs","test.log"),
  first_index=1,last_index=1,
  sim_repeat=2,
  stim_count=25,
  git_hash="UNKNOWN",
  params=joinpath(datadir,"params.jld2"),
  progressbar=false
)

