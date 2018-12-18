using StatsBase: fit, Histogram, weights
using Random

rms(x) = sqrt(mean(x.^2))
meansqr(x) = mean(x.^2)

function select_params(params;kwds...)
  condition = trues(size(params,1))
  for (var,val) in pairs(kwds)
    condition .&= ustrip.(abs.(params[var] .- val)) .<= 1e-2
  end
  params[condition,:pindex]
end

function select_data(df,params;kwds...)
  sel = select_params(params;kwds...)
  dfsel = df[in.(df.pindex,Ref(sel)),:]
  found = params[unique(dfsel.pindex),:]
  if size(found,1) < length(sel)
    @warn("Expected $(length(sel)) parameter entries. ",
          "\nInstead, only found entires: ",string(found),
          "\nKeyword selection: ",string(kwds))
  end
  dfsel[:st] = params.Δf[indexin(dfsel.pindex,params.pindex)]
  dfsel, params[sel,:]
end

function model_data(df,params;kwds...)
  df,params = select_data(df,params;kwds...)
  (stream=stream_summary(df,params), lengths=length_summary(df,params))
end

error_ratio(a,b=human_error()) = (a.stream / b.stream + a.lengths / b.lengths)/2
function model_error(df::DataFrame,params::DataFrame;kwds...)
  mdata = model_data(df,params;kwds...)
  model_error(mdata,HMEAN)
end

function model_error(data::NamedTuple,mean::NamedTuple)
  str_error = mean_bysid(data.stream) do str
    stream_rms(str,mean.stream)^2
  end |> sqrt
  len_error = ksstat(data.lengths,mean.lengths)

  (stream=str_error,lengths=len_error)
end

function mean_bysid(fn,data)
  total = 0.0
  n = 0
  for g in groupby(data,:sid)
    total += fn(g)
    n += 1
  end
  total / n
end

function stream_rms(data,mean)
  diffs = map(1:size(mean,1)) do row
    i = findfirst(isequal(mean.st[row]),data.st)
    if i isa Nothing
      missing
    else
      data[i,:streaming] - mean[row,:streaming]
    end
  end
  rms(coalesce.(diffs,0.0))
end

function ksstat(x,y)
  nx = length(x)
  ny = length(y)
  order = sortperm([x; y])
  diffs = cumsum(map(i -> i <= nx ? 1/nx : -1/ny,order))
  maximum(abs,diffs)
end

function findci(x;kwds...)
  if length(x) > 3
    if any(xi != x[1] for xi in x)
      cint = dbootconf(collect(skipmissing(x));kwds...)
      (mean=mean(skipmissing(x)), lowerc=cint[1], upperc=cint[2])
    else
      (mean=x[1], lowerc=x[1], upperc=x[1])
    end
  else
    (mean = mean(x), lowerc = minimum(x), upperc = maximum(x))
  end
end

function stream_summary(data,params;bound_threshold=0.8)
  result = by(data,[:st,:created]) do g
    DataFrame(streaming = streamprop(g.percepts,g.length,bound=true,
                                     threshold=bound_threshold))
  end
  if all(ismissing,result.streaming)
    result = by(data,[:st,:created]) do g
      DataFrame(streaming = streamprop(g.percepts,g.length,bound=false))
    end
  end
  result[:sid] = 0
  for g in groupby(result,:st)
    g[:sid] = 1:size(g,1)
  end
  sort!(result,(:sid,:st))

  result
end

function length_summary(data,params;Δf=6)
  pindex = params.pindex[params.Δf .== Δf]
  data[data.pindex .== pindex,:length]
end

function mean_streaming(df;findci=false)
  by(df,:st) do st
    if findci
      DataFrame([findci(st.streaming)])
    else
      DataFrame(streaming=mean(st.streaming))
    end
  end
end

function human_data(;resample=nothing)
  (stream=human_stream_data(),lengths=human_length_data(resample=resample))
end

function data_summarize(data)
  stream = by(data.stream,:st) do g
    DataFrame(streaming = mean(g.streaming))
  end
  lengths = data.lengths

  (stream=stream,lengths=lengths)
end

function human_error(;kwds...)
  str, len = human_error_by_sid(;kwds...)
  (stream=mean(str.x1), lengths=mean(len.x1))
end

function human_error_by_sid(;resample=1000,N=1)
  means,meanl = data_summarize(human_data())
  stream = human_stream_data()
  lengths = human_length_data(resample=resample,N=N)

  str_error = map(groupby(stream,:sid)) do str
    stream_rms(str,means)
  end |> combine
  len_error = map(groupby(lengths,:sid)) do len
    ksstat(len.lengths,meanl)
  end |> combine

  str_error, len_error
end

function human_stream_data()
  df = CSV.read(joinpath("..","analysis","context","stream_prop.csv"))
  rename!(df,:response => :streaming)
  sort!(df,(:sid,:st))
  df
end

const N_for_pressnitzer_hupe_2006 = 23
const pressnitzer_hupe_binsize = 1/6

function human_length_data(;resample=nothing,N = N_for_pressnitzer_hupe_2006)
  ph = CSV.read(joinpath("..","data","pressnitzer_hupe",
                         "pressnitzer_hupe_inferred.csv"))
  ph.length .+= pressnitzer_hupe_binsize.*(-0.5.+rand(size(ph,1)))
  lengths = collect(skipmissing(ph.length))

  if resample isa Nothing
    lengths
  else
    nsamples = div(length(lengths),N)
    inds = dbootinds(lengths,numobsperresample=nsamples,numresample=resample)
    dfs = map(enumerate(inds)) do (i,inds)
      DataFrame(lengths = lengths[inds], sid = i)
    end

    vcat(dfs...)
  end
end

# THOUGHTS: in P&H 2006, the mean normalization is on a per-individual basis.
# This might incline us to use a mean on a per-simulation basis, but I think it
# makes more sense to do across all runs, because the simulation represents
# repeated measurements from the same "individual"
function normlength(x)
  x = log.(x)
  x ./= mean(x)
  s = std(x.-1)
  if !iszero(s)
    x ./= s
  end
  exp.(x)
end


function handlebound(fn,seconds;bound=true,threshold=0.8)
    if length(seconds) < 2
      fn(1:length(seconds))
    end

    if bound && (sum(seconds[2:end]) > threshold*sum(seconds))
      fn(2:length(seconds))
    else
      fn(1:length(seconds))
    end
end

function buildup_mean(buildup_data;delta,length)
  buildup = DataFrame(time=range(0,stop=length,step=delta))
  buildup[:value] = 0.0
  runs = groupby(buildup_data,:run)
  for run in runs
    j = 1
    ts = cumsum(run.length)
    for (i,t) in enumerate(buildup.time)
      while j <= Base.length(ts) && t > ts[j]
        j += 1
      end
      j <= Base.length(ts) || break
      buildup.value[i] += (buildup_data.response[j]-1)
    end
  end
  buildup.value ./= Base.length(runs)
  buildup
end

function streamprop(percepts,seconds;kwds...)
  handlebound(seconds;kwds...) do range
    sum(seconds[range][percepts[range] .== 2]) / sum(seconds[range])
  end
end

function stim_per_second(seconds;kwds...)
  handlebound(seconds;kwds...) do range
    length(range) / sum(seconds[range])
  end
end

const HMEAN = data_summarize(human_data())
const HERR = human_error()
