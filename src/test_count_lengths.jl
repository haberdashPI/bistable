using Pkg; Pkg.activate("..")
include("count_lengths.jl")

# TODO: why are 3st and 12st analyses so similar to one another in the tracking
# model?? it looks like their similarities are very close, do we need a notion
# of neighbor hood, or can we leverage the multiescale representation to manage
# this

dir = joinpath("..","data","count_lengths")
isdir(dir) || mkdir(dir)

# select a specific set of parameters
paramfile = joinpath("..","data","count_lengths","run_2018-09-12",
                     "individual_extremes_params.feather")
params = getparams(paramfile) do i,row
  abs(row[:t_c_m] - 32) < 1 &&
  abs(row[:t_c_a] - 5) < 1 &&
  row[:Δf] == 12
end

params[:Δf] = 12
params[:t_W_m_σ_t] = 7.0
params[:t_W_m_σ_ϕ] = 7.0
params[:t_W_m_c] = 4.0
result12 = bistable_model(200, params, "settings.toml", interactive=true,
                          progressbar=false)
params[:Δf] = 6
result6 = bistable_model(200, params, "settings.toml", interactive=true,
                         progressbar=false)
params[:Δf] = 3
result3 = bistable_model(200, params, "settings.toml", interactive=true,
                         progressbar=false)
params[:Δf] = 1
result1 = bistable_model(200, params, "settings.toml", interactive=true,
                         progressbar=false)
# using RCall
alert()

# current plan: develop some intuition abou the role the scales

# are playing in evaluating the neighborhood,
# then think about what changes might help

# theory: the larger scales favor grouping (because the responses appear closer)
# test the theory by eliminating the smaller scales

# that seems to be about right...  changing the scale parameters seems to move
# us closer to expected behavior (suggests the underlying issue is the weighting
# and/or normalization of the scales)
#
# but it's not perfect, there's something non-monotonic about it...
#
# plan:
#
# I want to see all of the steps to get to the probability
#
# maybe start with idealized frames and their distances
# move to tracking of those idealized frames

# also need to consider that the greedy grouping isn't
# working well

# this is always run separately from the above code
# and they should probably be split into two files
# call bistable model here
datadir = joinpath(data_dir,"count_lengths","run_2018-09-12")
count_lengths(
  datadir=datadir,
  first_index=1,last_index=150,
  sim_repeat=10,
  stim_count=100,
  params=joinpath(datadir,"individual_extremes_params.feather"),
  progressbar=false
)

