using ArgParse
include(joinpath(@__DIR__,"count_lengths.jl"))

parse_settings = ArgParseSettings()
@add_arg_table parse_settings begin
  "first_index"
    help = "First parameter index to test."
    required = false
    arg_type = Int
    default = 1
  "last_index"
    help = "Last parameter index to test."
    required = false
    arg_type = Int
    default = 0
  "--params"
    help = "The file specifying a databse of parameters to explore."
    requried = false
    arg_type = String
    default = joinpath(@__DIR__,"params.jld2")
  "--repeat", "-r"
    help = "How many times to repeat simulation for each parameter."
    required = false
    arg_type = Int
    default = 1
  "--stim_count", "-c"
    help = "How many counts of the ab stimulus to present in each simulation."
    required = false
    arg_type = Int
    default = 50
  "--datadir", "-d"
    help = "The directory to save results to"
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"..","..","..","data","count_lengths")
  "--logfile", "-l"
    help = "The file to log information to"
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"..","..","..","data","count_lengths","run.log")
  "--scale_start"
    help = "Power of 2 to start cortical scales at."
    required = false
    arg_type = Float64
    default = -1.0
  "--scale_stop"
    help = "Power of 2 to stop cortical scales at."
    required = false
    arg_type = Float64
    default = 2.0
  "--scale_N"
    help = "Number of cortical scales to use."
    required = false
    arg_type = Int
    default = 12
  "--rate_start"
    help = "Power of 2 to start cortical rates at."
    required = false
    arg_type = Float64
    default = 1.0
  "--rate_stop"
    help = "Power of 2 to end cortical rates at."
    required = false
    arg_type = Float64
    default = 5.0
  "--rate_N"
    help = "Number of cortical rates to use."
    required = false
    arg_type = Int
    default = 4
  "--min_Hz"
    help = "Minimum Hz to model in the rates."
    required = false
    arg_type = Float64
    default = 400.0
  "--max_Hz"
    help = "Maximum Hz to model in the rates."
    required = false
    arg_type = Float64
    default = 1200.0
  "--nmf_K"
    help = "Number of NMF components to compute."
    required = false
    arg_type = Int
    default = 3
  "--nmf_window"
    help = "NMF analysis window size (in ms)."
    required = false
    arg_type = Float64
    default = 200.0
  "--nmf_delta"
    help = "NMF analysis window step (in ms)."
    required = false
    arg_type = Float64
    default = 100.0
  "--nmf_itr"
    help = "Maximum number of iterations of NMF analysis per window."
    required = false
    arg_type = Int
    default = 50
  "--track_prior_freq_N"
    help = "N for prior of novel source frequency. The higher N the stronger the"*
           " priors effect."
    required = false
    arg_type = Float64
    default = 2.0
  "--track_prior_freq_bias"
    help = "The bias for the prior of novel source frequency. bias = 0 no new"*
           "sources, bias 1 = all new sources."
    required = false
    arg_type = Float64
    default = 0.1
  "--track_max_sources"
    help = "The maximum number of sources to track."
    required = false
    arg_type = Int
    default = 3
  "--track_tc_start"
    help = "The first source tracking time constant."
    required = false
    arg_type = Float64
    default = 0.1
  "--track_tc_end"
    help = "The last source tracking time constant."
    required = false
    arg_type = Float64
    default = 0.25
  "--track_tc_N"
    help = "The number of source tracking time constants."
    required = false
    arg_type = Int
    default = 2
  "--track_zero_thresh"
    help = "The threshold for treating responses as zero during source"*
    "tracking."
    required = false
    arg_type = Float64
    default = 1e-1
  "--track_early_C"
    help = "The amount of source data in seconds to use as a prior."
    required = false
    arg_type = Float64
    default = 2.0
  "--track_prior_strength"
    help = "The strength of the source prior. The larger the more influence on"*
    "source models the prior has."
    required = false
    arg_type = Float64
    default = 10.0
end
args = parse_args(parse_settings)

count_lengths(args)
