include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

# NOTE: for freq W_m_σ is 2.0
# for scales W_m_σ is 16.0
params = Dict(
  :delta_t    => 240ms,                         :delta_f   => 3,
  :standard_f => 500Hz,                         :condition => :freqs,
  :c_x        => 3.0,                           :τ_x       => 500ms,
  :c_σ        => 0.2,                           :τ_σ       => 500ms,
  :c_a        => 6,                             :τ_a       => 3s,
  :c_m        => 65,                            :τ_m       => 350ms,
  :W_m_σ      => 5.6,                           :W_m_c     => 6.0,
 )

settings = TOML.parsefile("settings_2018-07-02.toml")
# settings = TOML.parsefile("fast_settings.toml")
# result = bistable_model(20, params, settings, interactive=true)
# alert()

Logging.configure(level=INFO)
count_lengths(
  1,2,
  stim_count=2,
  params="params_2018-08-04.feather",
  git_hash="UNKNOWN",
  settingsfile="fast_settings.toml",
  progressbar=false
)

