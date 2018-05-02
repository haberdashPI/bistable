using Query
using FileIO
using Feather
using CSV
using DataFrames
using ProgressMeter
using DataFrames

push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence

include("util/stim.jl")
include("util/peaks.jl")
include("util/lengths.jl")
include("util/biscales.jl")
include("util/threshold.jl")

function percept_lengths(counts,minlen=0.5s)
  thresh = ustrip(minlen/Δt(counts))
  lens,vals = findlengths(counts)
  slens = lens * ustrip(Δt(counts))

  mergelengths(slens,vals,ustrip(minlen))
end

# TODO:
# - make sure there aren't any aggregious type instabilities
#   in the code that is running most frequently
# - use Logging to report status of program

function count_lengths_helper(x,methods,params)
  y = bistable_scales(x,params)
  vals = map(methods) do name_method
    name,method = name_method
    try
      len,stim = y |> method |> percept_lengths
      name => DataFrame(length = len,stimulus = stim,method = string(name))
    catch e
      if e isa ResponseOverflow
        info("Skipping a simulation for $prefix due to parameter overflow",
             " (maximum value = $(e.val))")
        name => DataFrame()
      else
        rethrow(e)
      end
    end
  end

  vcat(values(vals)...)
end

function count_lengths(prefix,dir,N,M,methods,params)
  x = ab(120ms,120ms,1,50,500Hz,6) |> normpower |> amplify(-10dB)

  for i = 1:N
    rows = (count_lengths_helper(x,methods,params) for i in 1:M)
    name = @sprintf("counts_length_%s_thread%03d_%s.feather",
                    prefix,i,string(Date(now())))
    df = vcat(rows...)
    Feather.write(joinpath(dir,name),df)
  end
end
