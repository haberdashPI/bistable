using TOML
using Pkg; Pkg.activate("..")
include("count_lengths.jl")

# select a specific set of parameters
paramfile = joinpath("..","data","count_lengths","run_2018-09-12",
                     "individual_extremes_params.feather")
params = getparams(paramfile) do i,row
  abs(row[:t_c_m] - 32) < 1 &&
  abs(row[:t_c_a] - 5) < 1 &&
  row[:Δf] == 12
end

params[:t_W_m_σ_t] = 7.0
params[:t_W_m_σ_ϕ] = 7.0
params[:t_W_m_c] = 4.0

settings = TOML.parsefile("settings.toml")

result = []
for d in [1,3,6,12]
  params[:Δf] = d
  push!(result,bistable_model(15,params,settings,progressbar=false,
                              intermediate_results=true))
end

# this makes little sense to me,
# 1. why are the split models more likely?
# 2. why are there fewer split models for the smaller
# Δf

# is this something weird about considering 3 sources?
settings["track"]["analyze"]["max_sources"] = 2

result = []
for d in [1,3,6,12]
  params[:Δf] = d
  push!(result,bistable_model(15,params,settings,progressbar=false,
                              intermediate_results=true))
end

# no

# next step: try out a few ways to normalize
# my last attempt did not work well,
# now I'm only going to normalize up to some limit (vectors that
# are too small don't get fully enlarged)
#
# if that fails, take a more careful look at the distances
