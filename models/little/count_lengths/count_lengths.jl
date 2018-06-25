using Logging
using FileIO
using DataFrames
using DataStructures
using Feather
using JLD2
using Parameters

push!(LOAD_PATH,joinpath(@__DIR__,"..","packages"))
using AuditoryModel
using AuditoryCoherence

include(joinpath(@__DIR__,"..","util","stim.jl"))
include(joinpath(@__DIR__,"..","util","peaks.jl"))
include(joinpath(@__DIR__,"..","util","lengths.jl"))
include(joinpath(@__DIR__,"..","util","biscales.jl"))
include(joinpath(@__DIR__,"..","util","threshold.jl"))

@with_kw struct CountLength
  length::Float64 = 0.0
  stimulus::Int = 0
  error_code::Int = 0
  pindex::Int
  created::DateTime = DateTime()
end

function for_count_lengths(fn,dir)
  for file in readdir(dir)
    if ismatch(r"jld2$",file)
      jldopen(joinpath(dir,file),"r") do stream
        for key in keys(stream)
          fn(stream[key])
        end
      end
    end
  end
end

const RESPONSE_OVERFLOW=1
function count_lengths_helper(x,param_index,start_time,params,resolution)
  try
    len,stim = bistable_model(x,params,args) |> percept_lengths
    map(zip(len,stim)) do len_stim
      len,stim = len_stim
      CountLength(length=len,stimulus=stim,pindex=param_index,created=start_time)
    end
    vcat(values(vals)...)
  catch e
    if e isa ResponseOverflow
      [CountLength(error_code=RESPONSE_OVERFLOW,created=start_time,
                   pindex=param_index)]
    else
      rethrow(e)
    end
  end
end

totime(x) = x.*ms
tofreq(x) = x.*Hz
function count_lengths_runner(args)
  dir = abspath(args["datadir"])
  isdir(dir) || mkdir(dir)

  info("Loading parameters from "*args["params"])
  params = Feather.read(args["params"],transforms = Dict(
    "τ_σ" => totime,
    "τ_m" => totime,
    "τ_a" => totime,
    "τ_x" => totime,
    "τ_n" => totime
    "condition" => x -> Symbol.(x)
  ))
  first_index = args["first_index"]
  if first_index > nrow(params)
    err("First parameter index ($first_index) is outside the range of",
        " parameters.")
  end
  last_index = clamp(args["last_index"],first_index,nrow(params))
  indices = first_index:last_index
  info("Reading parameters for indices $(indices)")

  scales = cycoct.*2.0.^linspace(args["scale_start"],args["scale_stop"],
                                 args["scale_N"])

  stim = ab(120ms,120ms,1,args["stim_count"],500Hz,6) |>
    normpower |> amplify(-10dB)
  stim_resp = cortical(audiospect(stim,progressbar=false),progressbar=false,
                       scales=scales)
  info("Generated stimulus with cortical scales from 2^$(args["scale_start"])",
       " to 2^$(args["scale_stop"]) in $(args["scale_N"]) steps.")

  @assert log10(nrow(params)) < 6

  info("Total threads: $(Threads.nthreads())")
  sim_repeat = args["repeat"]
  #=Threads.@threads=# for i in repeat(indices,inner=sim_repeat)
    rows = count_lengths_helper(stim_resp,methods,i,now(),
                                Dict(k => params[i,k] for k in names(params)),
                                args)

    name = @sprintf("results_params%06d_%06d_t%02d.jld2",
                    first_index,last_index,Threads.threadid())
    filename = joinpath(dir,name)
    jldopen(filename,"a+") do file
      count = length(keys(file))
      file[@sprintf("rows%02d",count)] = rows
    end

    if Threads.threadid() == 1
      info("Completed a run for paramter $i.")
      info("Saved $(length(rows)) to $name.")
    end
  end
  info("DONE")
  info("------------------------------------------------------------")
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
