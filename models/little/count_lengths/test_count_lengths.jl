include("count_lengths.jl")

dir = joinpath("..","..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

# NOTE: for freq W_m_σ is 2.0
# for scales W_m_σ is 16.0
params = Dict(
    :delta_t    => 240ms, :delta_f   => 3,
    :standard_f => 500Hz, :condition => :freqs,
    :τ_x        => 500ms, :c_x       => 3.0,
    :W_m_σ      => 2.0,   :W_m_c     => 6.0,
    :τ_m        => 350ms, :c_m       => 30,
    :τ_a        => 3s,    :c_a       => 10,
    :τ_σ        => 500ms, :c_σ       => 0.3
)

# ODD: delta_f for 2 now seems to be always splitting which is totally weird.
# Why is this? Is this specifically about the :track condition? or something
# 'minor' that I changed

# I do believe this is because of the :tracks vs. the :scales
# I should think a bit more about why this is, but should alos
# probably just run it, and see what happens

# DONE: test tracks for (12,0.5,6,3)

# NEXT: test (scale) with (3,0.5,12) and freq with (0.5,)3,12

# TODO: remove extra data saved in count_lenghts

settings = TOML.parsefile("settings_2018-07-02.toml")
# settings = TOML.parsefile("fast_settings.toml")
result = bistable_model(40, params, settings, interactive=true)
alert()

Logging.configure(level=INFO)
count_lengths(
  1,2,
  stim_count=10,
  git_hash="UNKNOWN",
  settingsfile="fast_settings.toml",
  progressbar=false
)

