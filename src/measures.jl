using StatsBase: fit, Histogram, weights
using Random

rms(x) = sqrt(mean(x.^2))
meansqr(x) = mean(x.^2)

function select_params(params;kwds...)
  condition = trues(size(params,1))
  for (var,val) in pairs(kwds)
    condition .&= abs.(params[var] .- val) .<= 0.1
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
  dfsel, params[sel,:]
end

function model_data(df,params;kwds...)
  df,params = select_data(df,params;kwds...)
  df[:st] = params.Δf[indexin(df.pindex,params.pindex)]
  (stream=stream_summary(df,params), lengths=length_summary(df,params))
end

const HMEAN = data_summarize(human_data())
const HERR = human_error()
error_ratio(a,b=HERR) = (a.stream / b.stream + a.lengths / b.lengths)/2
function model_error(df::DataFrame,params::DataFrame;kwds...)
  mdata = model_data(df,params;kwds...)
  model_error(mdata,HMEAN)
end

function model_error(data::NamedTuple,mean::NamedTuple)
  str_error = mean_bysid(data.stream) do str
    stream_rms(str,mean.stream)^2
  end |> sqrt
  len_error = length_rms(data.lengths,mean.lengths)

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

function length_rms(len,mean)
  dens = length_hist(len)
  rms(mean .- dens)
end

# manually selected range to ensure reasonable bin size given the bin size of
# the original data from P&H 2006
const normalized_hist_range = 0:0.33333:20
function normweights(x)
  total = sum(x)
  total > 0 ? x ./ total : x
end

function length_hist(len)
  len = collect(skipmissing(len))
  hist = fit(Histogram,len,weights(len),normalized_hist_range)
  normweights(hist.weights)
end

function findci(x;kwds...)
  if length(x) > 3
    if any(x != x[1] for x in x)
      cint = dbootconf(x;kwds...)
      (mean=mean(x), lowerc=cint[1], upperc=cint[2])
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
  lengths = length_hist(data.lengths)

  (stream=stream,lengths=lengths)
end

function human_error(;resample=1000)
  means,meanl = data_summarize(human_data())
  stream,lengths = human_data(resample=resample)

  str_error = mean_bysid(stream) do str
    stream_rms(str,means)^2
  end |> sqrt

  len_error = mean_bysid(lengths) do len
    length_rms(len.lengths,meanl)^2
  end |> sqrt

  (stream=str_error,lengths=len_error)
end

function human_stream_data()
  df = CSV.read(joinpath("..","analysis","context","stream_prop.csv"))
  rename!(df,:response => :streaming)
  sort!(df,(:sid,:st))
  df
end

const N_for_pressnitzer_hupe_2006 = 23
function human_length_data(;resample=nothing)
  ph = CSV.read(joinpath("..","data","pressnitzer_hupe",
                         "pressnitzer_hupe_inferred.csv"))
  if resample isa Nothing
    ph.length
  else
    lengths = collect(skipmissing(ph.length))
    nsamples = div(length(lengths),N_for_pressnitzer_hupe_2006)
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
    if length(seconds) < 3
      fn(1:length(seconds))
    end

    if bound && (sum(seconds[2:end-1]) > threshold*sum(seconds))
      fn(2:length(seconds)-1)
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
