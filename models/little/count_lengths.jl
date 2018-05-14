using Logging
using FileIO
using Feather
#  using CSV
using DataFrames
#  using ProgressMeter
using Feather
push!(LOAD_PATH,joinpath(@__DIR__,"packages"))
using AuditoryModel
using AuditoryCoherence

include(joinpath(@__DIR__,"util","stim.jl"))
include(joinpath(@__DIR__,"util","peaks.jl"))
include(joinpath(@__DIR__,"util","lengths.jl"))
include(joinpath(@__DIR__,"util","biscales.jl"))
include(joinpath(@__DIR__,"util","threshold.jl"))

# TODO:
# - make sure there aren't any aggregious type instabilities
#   in the code that is running most frequently
# - use Logging to report status of program

const RESPONSE_OVERFLOW=1
function count_lengths_helper(x,methods,params)
  try
    y = bistable_scales(x,params)
    vals = map(methods) do name_method
      name,method = name_method
      len,stim = y |> method |> percept_lengths
      name => DataFrame(length = len,stimulus = stim,method = string(name),
                        error_code=0)
    end
    vcat(values(vals)...)
  catch e
    if e isa ResponseOverflow
      DataFrame(length=0,stimulus=0,method="NONE",error_code=RESPONSE_OVERFLOW)
    else
      rethrow(e)
    end
  end
end

function count_lengths_runner(args)
  dir = abspath(args["datadir"])
  isdir(dir) || mkdir(dir)

  params = load(joinpath(@__DIR__,"params.jld2"),"df")
  first_index = args["first_index"]
  if first_index > nrow(params)
    err("First parameter index ($first_index) is outside the range of",
        " parameters.")
  end
  last_index = clamp(args["last_index"],first_index,nrow(params))
  indices = first_index:last_index
  info("Reading parameters for indices $(indices)")

  sim_repeat = args["repeat"]
  # TODO: change these counts back to higher values and make
  # use of threads (julia will also need to be run so that it
  # has more than one thread, check how many I can actually
  # use).

  methods = Dict(
    :threshold => x -> source_count_by_threshold(x,window=1s,delta=0.25s,
                                                 cutoff=2cycoct,buildup=1s)
  # :peaks => x -> source_count_by_peaks(x,window=1s,delta=0.25s,buildup=1s),
  )
  info("Testing with the method$(length(methods) > 1 ? "s" : "") "*
       join(map(string,keys(methods)),", "," and ")*".")

  scales = cycoct.*2.0.^linspace(args["scale_start"],args["scale_stop"],
                                 args["scale_N"])

  stim = ab(120ms,120ms,1,args["stim_count"],500Hz,6) |>
    normpower |> amplify(-10dB)
  stim_resp = cortical(audiospect(stim,progressbar=false),progressbar=false,
                       scales=scales)
  info("Generated stimulus with cortical scales from 2^$(args["scale_start"])",
       " to 2^$(args["scale_stop"]) in $(args["scale_N"]) steps.")

  results = Array{DataFrame}(length(indices))
  info("Total threads: $(Threads.nthreads())")
  Threads.@threads for i in indices
    rows = vcat((count_lengths_helper(stim_resp,methods,
                                      Dict(k => params[i,k] for k in
                                           names(params)))
                 for repeat in 1:sim_repeat)...)
    rows[:,:param_index] = i
    results[i] = rows
  end

  @assert log10(nrow(params)) < 6
  name = @sprintf("results_params%06d_%06d.feather",first_index,last_index)
  filename = joinpath(dir,name)
  Feather.write(filename,vcat(results...))
  info("Wrote results to $name")
end

function count_lengths(args)
  info("Results will be saved to $(args["datadir"]).")
  info("All logging data will be saved to $(args["logfile"]).")
  Logging.configure(output=open(args["logfile"], "a"),level=INFO)
  try
    count_lengths_runner(args)
  catch ex
    str = sprint(io->Base.show_backtrace(io, catch_backtrace()))
    err("$ex: $str")
  end
end
