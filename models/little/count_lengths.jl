using Logging
using FileIO
using DataFrames
using DataStructures

push!(LOAD_PATH,joinpath(@__DIR__,"packages"))
using AuditoryModel
using AuditoryCoherence

include(joinpath(@__DIR__,"util","stim.jl"))
include(joinpath(@__DIR__,"util","peaks.jl"))
include(joinpath(@__DIR__,"util","lengths.jl"))
include(joinpath(@__DIR__,"util","biscales.jl"))
include(joinpath(@__DIR__,"util","threshold.jl"))

struct CountLength
  length::Float64
  stimulus::UInt32
  method::UInt32
  error_code::UInt32
  param_index::UInt32
end

function Base.write(io::IO,row::CountLength)
  write(io,row.length)
  write(io,row.stimulus)
  write(io,row.method)
  write(io,row.error_code)
  write(io,row.param_index)
end
function saverows(file,rows::Array{CountLength})
  open(file,"a+") do io
    for row in rows; write(io,row); end
  end
end

function Base.read(io::IO,::Type{CountLength})
  CountLength(read(io,Float64),read(io,UInt32),read(io,UInt32),
              read(io,UInt32),read(io,UInt32))
end
function loadrows(file)
  rows = Array{CountLength}(0)
  open(file,"r") do io
    while !eof(io)
      push!(rows,Base.read(io,CountLength))
    end
  end
  rows
end

const RESPONSE_OVERFLOW=1
const method_index = Dict(:threshold => 1,:peaks => 2,:cohere => 3)
function count_lengths_helper(x,methods,param_index,params)
  try
    y = bistable_scales(x,params)
    vals = map(methods) do name_method
      name,method = name_method
      len,stim = y |> method |> percept_lengths
      name => CountLength.(len,stim,method_index[name],0,param_index)
    end
    vcat(values(vals)...)
  catch e
    if e isa ResponseOverflow
      [CountLength(0,0,0,RESPONSE_OVERFLOW,param_index)]
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

  methods = Dict(
    #  :threshold => x -> source_count_by_threshold(x,window=1s,delta=0.25s,
                                                 #  cutoff=2cycoct,buildup=1s)
    :peaks => x -> source_count_by_peaks(x,window=1s,delta=0.25s,buildup=1s)
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

  @assert log10(nrow(params)) < 6

  info("Total threads: $(Threads.nthreads())")
  sim_repeat = args["repeat"]
  Threads.@threads for i in repeat(indices,inner=sim_repeat)
    rows = count_lengths_helper(stim_resp,methods,i,
                                Dict(k => params[i,k] for k in names(params)))

    name = @sprintf("results_params%06d_%06d_t%02d.clbin",
                    first_index,last_index,Threads.threadid())
    filename = joinpath(dir,name)
    saverows(filename,rows)
    if Threads.threadid() == 1
      info("Completed a run for paramter $i.")
    end
  end
  info("DONE")
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
