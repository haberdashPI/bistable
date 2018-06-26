include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

params = Dict(
    :delta_t    => 240ms, :delta_f   => 6,
    :standard_f => 500Hz, :condition => :scales,
    :τ_x        => 500ms, :c_x       => 3.0,
    :W_m_σ      => 15.0,  :W_m_c     => 6.0,
    :τ_m        => 350ms, :c_m       => 30,
    :τ_a        => 3s,    :c_a       => 10,
    :τ_σ        => 500ms, :c_σ       => 0.3,
)

count_lengths(
  1,2,
  stim_count=10,
  git_hash="UNKNOWN",
  settingsfile="fast_settings.toml",
  progressbar=true
)

