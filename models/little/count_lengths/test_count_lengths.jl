include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

# NOTE: for freq W_m_σ is 2.0
# for scales W_m_σ is 16.0
params = Dict(
    :Δt         => 240ms, :Δf        => 3,
    :f          => 500Hz, :condition => :track,
    :τ_x        => 500ms, :c_x       => 3.0,
    :W_m_σ_t    => 6.0,   :W_m_σ_ϕ   => 8.0,
    :W_m_c      => 6.0,
    :τ_m        => 350ms, :c_m       => 100,
    :τ_a        => 3s,    :c_a       => 6,
    :τ_σ        => 500ms, :c_σ       => 0.2
)
# in progress, getting the new revisions to work
# CURRENT ISSUE: bistable track no longer bistable

# check the output to make sure it's reasonable
# see how much faster this version is

# CURRENT: bistability or tracks seems to be at least partially working, time to
# speed things up.

# problem - in this configuration the fused
# percept just isn't all that plausible
# (try looking at time constants again????)

# validate on all three stimuli and all three
# bistable conditions
settings = TOML.parsefile("settings_2018-08-18.toml")
# settings["percept_lengths"]["threshold"] = 1.6
# settings = TOML.parsefile("fast_settings.toml")
@time result = bistable_model(40, params, settings, interactive=true,
                              progressbar=false)
alert()

settings["percept_lengths"]["min_length_ms"] = 100
context = ab(120ms,120ms,1,40,500Hz,12)
x = ab(120ms,120ms,1,40,500Hz,3)
stim = [context; x] |> normpower |> amplify(-10dB)
context_result = bistable_model(stim,params,settings,interactive=true)

x = ab(120ms,120ms,1,40,500Hz,3) |> normpower |> amplify(-10dB)
test_alone_result = bistable_model(x,params,settings,interactive=true)
alert()

# Logging.configure(level=INFO)
# count_lengths(
#   1,2,
#   stim_count=2,
#   params="freq_params_2018-08-05.feather",
#   git_hash="UNKNOWN",
#   settingsfile="fast_settings.toml",
#   progressbar=false
# )

