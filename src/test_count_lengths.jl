using Pkg; Pkg.activate("..")
include("count_lengths.jl")

datadir = joinpath(data_dir,"count_lengths","run_2018-10-03")
count_lengths(
  datadir=datadir,
  first_index=1,last_index=1,
  sim_repeat=2,
  stim_count=5,
  params=joinpath(datadir,"individual_levels_params.feather"),
  progressbar=false
)

