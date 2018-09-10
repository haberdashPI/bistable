
include("count_lengths.jl")

# TODO: test how well analysis works when there is noise in each level

dir = joinpath("..","data","count_lengths")
isdir(dir) || mkdir(dir)

# NOTE: for freq W_m_σ is 2.0
# for scales W_m_σ is 16.0

params = Dict(
  :Δt         => 240ms, :Δf        => 3,
  :f          => 500Hz, :condition => :track,
  :τ_x        => 500ms, :c_x       => 3.0,
  :f_c_a => 0, :f_c_m => 0, :f_c_σ => 0,
  :s_c_a => 0, :s_c_m => 0, :s_c_σ => 0,
  :t_W_m_σ      => 15.0, #5.0
  :t_W_m_σ_t    => 8.0,   :t_W_m_σ_ϕ   => 6.0,
  :t_W_m_c      => 6.0,
  :t_τ_m        => 350ms, :t_c_m       => 100,
  :t_τ_a        => 3s,    :t_c_a       => 6,
  :t_τ_σ        => 500ms, :t_c_σ       => 0.2
 )

# call bistable model here

count_lengths(
  first_index=1,last_index=75,
  logfile="individual_simulation.log",
  sim_repeat=10,
  stim_count=100,
  params="individual_params_2018-09-09.feather",
  settings="settings_2018-09-07.toml",
  progressbar=false
)

