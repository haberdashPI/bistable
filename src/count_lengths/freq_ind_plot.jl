using Feather
using JLD2
using FileIO
include("count_lengths.jl")

dir = "../../../plots/freq_ind_$(Date(now()))"
isdir(dir) || mkdir(dir)

param_file = joinpath("..","..","..","data","count_lengths",
                      "params_2018-07-16.feather")

all_params = Feather.read(param_file,transforms = Dict{String,Function}(
    "τ_x" => x -> totime.(x),
    "τ_σ" => x -> totime.(x),
    "τ_a" => x -> totime.(x),
    "τ_m" => x -> totime.(x),
    "Δt" => x -> totime.(x),
    "f" => x -> tofreq.(x),
    "condition" => x -> Symbol.(x)
  ))

########################################
# what's going on with fusing 12 st percepts?

sindex = 3183
params = Dict(k => all_params[sindex,k] for k in names(all_params))
params[:θ] = 2.1
settings = TOML.parsefile("settings_2018-07-02.toml")

result = bistable_model(40,params,settings,interactive=true)

p = rplot(result[4])
R"""
ggsave($(joinpath(dir,"freq_12st_m27_a10_sig5.6.pdf")),$p)
"""

########################################
# what's going on with the same parameters at 3st?
# basically, there's less competition at the closer response values

sindex = 2967
params = Dict(k => all_params[sindex,k] for k in names(all_params))
params[:θ] = 2.1
settings = TOML.parsefile("settings_2018-07-02.toml")

result = bistable_model(40,params,settings,interactive=true)

p = rplot(result[4])
R"""
ggsave($(joinpath(dir,"freq_3st_m27_a10_sig5.6.pdf")),$p)
"""

########################################
# what's the best bistability look like for freq?

sindex = 2444
params = Dict(k => all_params[sindex,k] for k in names(all_params))
params[:θ] = 2.1
settings = TOML.parsefile("settings_2018-07-02.toml")

result = bistable_model(40,params,settings,interactive=true)

p = rplot(result[4])
R"""
ggsave($(joinpath(dir,"freq_6st_selective.pdf")),$p)
"""

########################################
# let's look at bistability in a place that might be working
# by our newer definition

sindex = 4256
params = Dict(k => all_params[sindex,k] for k in names(all_params))
params[:θ] = 2.1
settings = TOML.parsefile("settings_2018-07-02.toml")

cache = "cached_model.jld2"
if isfile(cache)
  result = load(cache,"result")
else
  result = bistable_model(40,params,settings,interactive=true)
  save(cache,"result",result)
end

crmask,norm = mask(result[7],result[3],settings,progressbar=true)

# now with 12 st
sindex = 4472
params = Dict(k => all_params[sindex,k] for k in names(all_params))
params[:θ] = 2.1
settings = TOML.parsefile("settings_2018-07-02.toml")

cache = "cached_model_12st.jld2"
if isfile(cache)
  result = load(cache,"result")
else
  result = bistable_model(40,params,settings,interactive=true)
  save(cache,"result",result)
end

crmask,norm = mask(result[7],result[3],settings,progressbar=true)
crmask2,norm = mask(result[7],result[3],settings,progressbar=true,order=2)


# now with 0.5 st
sindex = 4040
params = Dict(k => all_params[sindex,k] for k in names(all_params))
params[:θ] = 2.1
settings = TOML.parsefile("settings_2018-07-02.toml")

cache = "cached_model_0.5st.jld2"
if isfile(cache)
  result = load(cache,"result")
else
  result = bistable_model(40,params,settings,interactive=true)
  save(cache,"result",result)
end

# TODO: create some graphs for mounya???
