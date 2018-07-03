using Logging
using FileIO
using DataFrames
using DataStructures
using Feather
using JLD2
using Parameters
using TOML

push!(LOAD_PATH,joinpath(@__DIR__,"..","packages"))
using AuditoryModel
using AuditoryCoherence

include(joinpath(@__DIR__,"..","util","stim.jl"))
include(joinpath(@__DIR__,"..","util","peaks.jl"))
include(joinpath(@__DIR__,"..","util","lengths.jl"))
include(joinpath(@__DIR__,"..","util","bimodel.jl"))
include(joinpath(@__DIR__,"..","util","biapply.jl"))
include(joinpath(@__DIR__,"..","util","threshold.jl"))

@with_kw struct CountLength
  length::Float64
  stimulus::Int
  pindex::Int
  created::DateTime
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

function count_lengths(args::Dict)
  info("Results will be saved to $(args["datadir"]).")
  info("All logging data will be saved to $(args["logfile"]).")
  Logging.configure(output=open(args["logfile"], "a"),level=INFO)
  try
    count_lengths(args["first_index"],args["last_index"],
                  params=args["params"],
                  git_hash=args["git_hash"],
                  sim_repeat=args["sim_repeat"],
                  stim_count=args["stim_count"],
                  datadir=args["datadir"],
                  logfile=args["logfile"],
                  settingsfile=args["settings"])
  catch ex
    str = sprint(io->Base.show_backtrace(io, catch_backtrace()))
    err("$ex: $str")
  end
end

function read_git_hash()
  olddir = pwd()
  cd(@__DIR__)
  hash = readstring(`git rev-parse HEAD`)
  cd(olddir)

  hash
end

totime(x) = x.*ms
tofreq(x) = x.*Hz
const data_dir = joinpath(@__DIR__,"..","..","..","data")
function count_lengths(first_index,last_index;
                       params=joinpath(@__DIR__,"params.jld2"),
                       git_hash="DETECT",
                       sim_repeat=2,
                       stim_count=25,
                       datadir=joinpath(data_dir,"count_lengths"),
                       logfile=joinpath(data_dir,"count_lengths","run.log"),
                       settingsfile=joinpath(@__DIR__,"settings.toml"),
                       progressbar=false)
  dir = abspath(datadir)
  isdir(dir) || mkdir(dir)

  info("Source code hash: "*(git_hash == "DETECT" ? read_git_hash() : git_hash))

  info("Loading parameters from "*params)
  params = load(params,"params")
  first_index = first_index
  if first_index > nrow(params)
    err("First parameter index ($first_index) is outside the range of",
        " parameters.")
  end
  last_index = clamp(last_index,first_index,nrow(params))
  indices = first_index:last_index
  info("Reading parameters for indices $indices")

  settings = TOML.parsefile(settingsfile)
  info("Reading settings from file $settingsfile")

  @assert log10(nrow(params)) < 4

  for i in repeat(indices,inner=sim_repeat)
    start_time = now()
    params_dict = Dict(k => params[i,k] for k in names(params))
    len,stim = bistable_model(stim_count,params_dict,settings,
                              progressbar=progressbar)
    rows = map(zip(len,stim)) do len_stim
      len,stim = len_stim
      CountLength(length=len,stimulus=stim,pindex=i,created=start_time)
    end

    name = @sprintf("results_params%04d_%04d.jld2",first_index,last_index)
    filename = joinpath(dir,name)
    jldopen(filename,"a+") do file
      count = length(keys(file))
      file[@sprintf("rows%02d",count)] = rows
    end

    info("Completed a run for paramter $i.")
    info("Saved $(length(rows)) row(s) to $name.")
  end
  info("DONE")
  info("------------------------------------------------------------")
end

