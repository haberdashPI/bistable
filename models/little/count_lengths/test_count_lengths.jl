include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

params = Dict(
    :delta_t    => 240ms, :delta_f   => 6,
    :standard_f => 500Hz, :condition => :scales,
    :τ_x        => 500ms, :c_x       => 3.0,
    :W_m_σ      => 15.0,  :W_m_c     => 6.0,
    :τ_m        => 350ms, :c_m       => 30,
    :τ_a        => 3s,   :c_a       => 10,
    :τ_σ        => 500ms, :c_σ       => 0.3
)

freq_params = Dict(
    :delta_t    => 240ms, :delta_f   => 6,
    :standard_f => 500Hz, :condition => :freqs,
    :τ_x        => 500ms, :c_x       => 3.0,
    :W_m_σ      => 10.0,  :W_m_c     => 6.0,
    :τ_m        => 350ms, :c_m       => 42,
    :τ_a        => 3s,    :c_a       => 20,
    :τ_σ        => 500ms, :c_σ       => 0.2
)


# TODO: have a min denominator for the ratio
# to avoid the wacky thing happening with the freq params
# (then fix what's wrong with freq params, where it's just blanking out the
# signal)

settings = TOML.parsefile("settings_2018-07-02.toml")
# settings = TOML.parsefile("fast_settings.toml")
result = bistable_model(60,freq_params, settings, interactive=true)

Logging.configure(level=INFO)
count_lengths(
  1,2,
  stim_count=10,
  git_hash="UNKNOWN",
  settingsfile="fast_settings.toml",
  progressbar=false
)
alert()
