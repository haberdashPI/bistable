include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

# NOTE: for freq W_m_σ is 2.0
# for scales W_m_σ is 16.0
params = Dict(
  :delta_t    => 240ms,                         :delta_f   => 3,
  :standard_f => 500Hz,                         :condition => :scales_track,
  :s_c_x      => 3.0,                           :s_τ_x     => 500ms,
  :s_c_σ      => 0.2,                           :s_τ_σ     => 500ms,
  :s_c_a      => 10,                            :s_τ_a     => 3s,
  :s_c_m      => 27,                            :s_τ_m     => 350ms,
  :s_W_m_σ    => 15.0,                          :s_W_m_c   => 6.0,
  :t_c_x      => 3.0,                           :t_τ_x     => 500ms,
  :t_c_σ      => 0.2,                           :t_τ_σ     => 500ms,
  :t_c_a      => 10,                            :t_τ_a     => 3s,
  :t_c_m      => 27,                            :t_τ_m     => 350ms,
  :t_W_m_σ    => 15.0,                          :t_W_m_c   => 6.0,
  :θ          => 2.0
 )

# settings = TOML.parsefile("settings_2018-07-02.toml")
settings = TOML.parsefile("fast_settings.toml")
result = bistable_model(2, params, settings, interactive=true)
alert()

Logging.configure(level=INFO)
count_lengths(
  1,2,
  stim_count=2,
  git_hash="UNKNOWN",
  settingsfile="fast_settings.toml",
  progressbar=false
)

