using Logging
using AuditoryBistabilityLE
using ShammaModel
using Dates
using Printf
using JLD2
using Feather
using InteractiveUtils
using Unitful

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
      @info "Loading $file"
      jldopen(joinpath(dir,file),"r") do stream
        for param in keys(stream)
          if occursin(r"param[0-9]{2}",param)
            for run in keys(stream[param])
              entry = stream[param][run]
              fn(CountLength(entry["ratio"], entry["bratio"],
                             entry["pindex"], entry["created"]))
            end
          end
        end
      end
    end
  end
end

count_lengths(args::Dict) = count_lengths(;(Symbol(k) => v for (k,v) in args)...)

function count_lengths(;first_index,last_index,
                       datadir=joinpath(data_dir,"count_lengths"),
                       logfile=joinpath(datadir,"sim.log"),
                       params=joinpath(datadir,"params.feather"),
                       dataprefix="results",
                       git_hash="DETECT",
                       sim_repeat=2,
                       stim_count=25,
                       settings=joinpath(@__DIR__,"settings.toml"),
                       progressbar=false)
  setup_logging(logfile) do
    @info("Results will be saved to $datadir")
    count_lengths(first_index,last_index,logfile,datadir,dataprefix,
                  params,git_hash,sim_repeat,stim_count,
                  settings,progressbar)
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
const data_dir = if occursin("Mycroft",gethostname())
  joinpath(@__DIR__,"..","data")
else
  joinpath(homedir(),"work","dlittle","bistable_individual")
end

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

function getparams(filterfn,file)
  asdict(params) = Dict(k => params[1,k] for k in names(params))

  params = handle_units!(Feather.read(file))
  rows = filter(ri -> filterfn(ri,asdict(params[ri,:])),Base.axes(params,1))
  if length(rows) > 1
    @warn("Specification ambiguous, multiple rows match.")
  elseif length(rows) < 1
    error("No rows matching specification")
  else
    asdict(params[rows[1],:])
  end
end

function count_lengths(first_index,last_index,logfile,datadir,dataprefix,
                       params,git_hash,sim_repeat,stim_count,
                       settings,progressbar)
  dir = abspath(datadir)
  isdir(dir) || mkdir(dir)

  verbuf = IOBuffer()
  versioninfo(verbuf)
  @info String(take!(verbuf))
  @info "Source code hash: "*(git_hash == "DETECT" ? read_git_hash() : git_hash)

  @info "Loading parameters from "*params
  params = handle_units!(Feather.read(params))
  first_index = first_index
  if first_index > size(params,1)
    error("First parameter index ($first_index) is outside the range of",
          " parameters.")
  end
  last_index = clamp(last_index,first_index,size(params,1))
  indices = first_index:last_index
  @info "Reading parameters for indices $indices"
  @info "Reading settings from file $settings"

  @assert log10(size(params,2)) < 5
  name = @sprintf("%s_params%05d_%05d.jld2",dataprefix,first_index,
                  last_index)
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
      result = bistable_model(stim_count,params_dict,settings,
                              progressbar=progressbar)

      jldopen(filename,"a+") do file
        # if it hasn't yet been done, record the times at which the ratios are
        # computed (their sample rate differs, so we need to interpolate them
        # later on).
        if !haskey(file,"btimes_s")
          file["btimes_s"] = ustrip.(uconvert.(s,times(result.percepts.ratio)))
        end

        entry = @sprintf("param%05d",i)
        count = haskey(file,entry) ? length(keys(file[entry])) : 0
        file[@sprintf("param%05d/run%03d/lengths",i,count)] =
          Array(result.percepts.counts[1])
        file[@sprintf("param%05d/run%03d/percepts",i,count)] =
          Array(result.percepts.counts[2])
        file[@sprintf("param%05d/run%03d/mask",i,count)] =
          compress!(result.primary_source)
        file[@sprintf("param%05d/run%03d/pindex",i,count)] = i
        file[@sprintf("param%05d/run%03d/created",i,count)] = start_time

        @info "Completed simulation run $(count+1) for parameter $i."
      end

      @info "Run yielded ~$(length(result.percepts.counts[1])) percepts."
    end
  end
  @info "DONE"
end

