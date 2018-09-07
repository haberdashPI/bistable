using Logging
using AuditoryBistabilityLE
using ShammaModel
using Dates
using Printf
using JLD2
using Feather
using InteractiveUtils
include("logger.jl")

struct CountLength
  ratio::Array{Float64}
  bratio::Array{Float64}
  pindex::Int
  created::DateTime
end

function for_count_lengths(fn,dir)
  for file in readdir(dir)
    if ismatch(r"jld2$",file)
      jldopen(joinpath(dir,file),"r") do stream
        for key in keys(stream)
          fn(CountLength(stream[key]["ratio"],stream[key]["bratio"],
                         stream[key]["pindex"],stream[key]["created"]))
        end
      end
    end
  end
end

function count_lengths(args::Dict)
  @info "All output will be saved to '$(args["logfile"])'."
  open(args["logfile"],"a") do stream
    with_logger(DatedLogger(stream)) do
      redirect_stdout(stream) do
        redirect_stderr(stream) do
          try
            @info("Results will be saved to $(args["datadir"]).")
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
            @error("$ex: $str")
          end
          @info("------------------------------------------------------------")
        end
      end
    end
  end
end

function read_git_hash()
  olddir = pwd()
  cd(@__DIR__)
  hash = read(`git rev-parse HEAD`,String)
  cd(olddir)

  hash
end

totime(x) = x.*ms
tofreq(x) = x.*Hz
const data_dir = joinpath(@__DIR__,"..","data")

from_unit_table = Dict(r"τ" => totime, r"Δt" => totime, r"^f$" => tofreq)
function handle_units!(df)
  for col in names(df)
    for (pattern,fn) in pairs(from_unit_table)
      if occursin(pattern,string(col))
        df[col] = fn.(df[col])
      end
    end
  end
  df
end

function count_lengths(first_index,last_index;
                       params=joinpath(@__DIR__,"params.feather"),
                       git_hash="DETECT",
                       sim_repeat=2,
                       stim_count=25,
                       datadir=joinpath(data_dir,"count_lengths"),
                       logfile=joinpath(data_dir,"count_lengths","run.log"),
                       settingsfile=joinpath(@__DIR__,"settings.toml"),
                       progressbar=false)
  dir = abspath(datadir)
  isdir(dir) || mkdir(dir)

  verbuf = IOBuffer()
  versioninfo(verbuf)
  @info String(take!(verbuf))
  @info "Source code hash: "*(git_hash == "DETECT" ? read_git_hash() : git_hash)

  @info "Loading parameters from "*params
  params = handle_units!(Feather.read(params))
  first_index = first_index
  if first_index > size(params,2)
    err("First parameter index ($first_index) is outside the range of",
        " parameters.")
  end
  last_index = clamp(last_index,first_index,size(params,1))
  indices = first_index:last_index
  @info "Reading parameters for indices $indices"
  @info "Reading settings from file $settingsfile"

  @assert log10(size(params,2)) < 5
  name = @sprintf("results_params%05d_%05d.jld2",first_index,last_index)
  filename = joinpath(dir,name)

  for i in indices
    @info "Running simulations for parameter $i"
    # determine the number of simulations to run, accounting for any previously
    # run simulations
    # NOTE: assumes only a single process ever accesses this file
    num_repeats = if isfile(filename)
      jldopen(filename,"r") do file
        p = @sprintf("param%05d",i)
        if haskey(file,p)
          max(0,sim_repeat - length(keys(file[p])))
        else
          sim_repeat
        end
      end
    else
      sim_repeat
    end

    @info "$(sim_repeat - num_repeats) simulations run previosuly."
    @info "Running $(num_repeats) more simulations."

    for repeat in 1:num_repeats
      start_time = now()
      params_dict = Dict(k => params[i,k] for k in names(params))
      result = bistable_model(stim_count,params_dict,settingsfile,
                              progressbar=progressbar)

      jldopen(filename,"a+") do file
        count = length(keys(file))
        file[@sprintf("param%05d/run%03d/ratio",i,count)] =
          Array(result.percepts.sratio)
        file[@sprintf("param%05d/run%03d/bratio",i,count)] =
          Array(result.percepts.bratio)
        file[@sprintf("param%05d/run%03d/pindex",i,count)] = i
        file[@sprintf("param%05d/run%03d/created",i,count)] = start_time
      end

      @info "Completed a simulation for paramter $i."
      @info "Run yielded ~$(length(result.percepts.counts)) percepts."
    end
  end
  @info "DONE"
end

