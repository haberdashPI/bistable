using Pkg; Pkg.activate("..")
include("count_lengths.jl")

datadir = joinpath(data_dir,"count_lengths","run_2018-10-03")
count_lengths(
  datadir=datadir,
  first_index=1,last_index=150,
  sim_repeat=10,
  stim_count=100,
  params=joinpath(datadir,"individual_levels_params.feather"),
  progressbar=false
)

