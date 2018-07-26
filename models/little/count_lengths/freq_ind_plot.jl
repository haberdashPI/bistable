using Feather
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
    "delta_t" => x -> totime.(x),
    "standard_f" => x -> tofreq.(x),
    "condition" => x -> Symbol.(x)
  ))

########################################
# what's going on with the fusing 12 st percepts?

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

# not super compelling, is this just a bad run?
# (if so we probably need a seed)
# might have to check on subsequent stages of the pipeline
