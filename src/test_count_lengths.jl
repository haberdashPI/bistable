using Pkg; Pkg.activate(joinpath(@__DIR__,".."))
include("count_lengths.jl")

datadir = joinpath(data_dir,"run_2018-10-04")
count_lengths(
  datadir=datadir,
  first_index=1,last_index=1,
  sim_repeat=2,
  stim_count=5,
  git_hash="UNKNOWN",
  params=joinpath(datadir,"individual_levels_params.feather"),
  progressbar=false
)

