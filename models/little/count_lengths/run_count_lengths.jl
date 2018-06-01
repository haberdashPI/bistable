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
    help = "Second parameter index to test."
    required = false
    arg_type = Int
    default = 0
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
    default = joinpath(@__DIR__,"..","..","data","count_lengths")
  "--logfile", "-l"
    help = "The file to log information to"
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"..","..","data","count_lengths","run.log")
end
args = parse_args(parse_settings)

count_lengths(args)
