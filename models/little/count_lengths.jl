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

function percept_lengths(method,csa,minlen=1s)
  counts = method(csa)
  thresh = ustrip(minlen/Δt(counts))
  lens,vals = findlengths(counts)
  slens = lens * ustrip(Δt(counts))

  mergelengths(slens,vals,ustrip(minlen))
end

# TODO:
# - walk through new source_count_by_threshold
#   function and make sure it is doing the right thing
# - check through code to make sure there aren't any additional 'hidden'
#   parameters
# - make sure there aren't any aggregious type instabilities
#   in the code that is running most frequently
# - use Logging to report status of program

function count_lengths_helper(x,methods,params)
  map(methods) do name_method
    y = bistable_scales(x,params)
    name,method = name_method
    try
      len,stim = percept_lengths(method,y)
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
end

function count_lengths(prefix,dir,N,methods,params)
  x = ab(120ms,120ms,1,50,500Hz,6) |> normpower |> amplify(-10dB)

  for i = 1:N
    rows = count_lengths_helper(x,methods,params)
    name = @sprintf("counts_length_%s_thread%03d_%s.feather",
                    prefix,i,string(Date(now())))

    Feather.write(joinpath(dir,name),vcat(values(rows)...))
  end
end
