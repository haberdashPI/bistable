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

function setup_human_data()
  (stream=stream_dfh(),lengths=length_hdata())
end

function model_rms(df,params,dfh;return_parts=false,N=1000,kwds...)
  (streaming,lengths) =
  (stream_rms(df,params,dfh.stream,return_parts=true,N=N;kwds...),
   length_rms(df,params,dfh.lengths,return_parts=true;kwds...))

  if return_parts
    (rms = mean((streaming.rms,lengths.rms)),
     stream_rms = streaming.rms,
     stream_sd = streaming.sd,
     stream_mratio = streaming.mratio,
     length_rms = lengths.rms,
     length_sd = lengths.sd)
  else
    mean((streaming.rms,lengths.rms))
  end
end

function stream_df(df,params;findci=true,bound=true,bound_threshold=0.8,kwds...)
  dfs = stream_dfs(df,params;findci=findci,bound=bound,
                   bound_threshold=bound_threshold,kwds...)
  dfh = stream_dfh(findci=findci)
  delete!(dfh,:rms)
  vcat(dfs,dfh)
end

function stream_dfs(df,params;resample_N=missing,resampling=Colon(),
                    keep_simulations=false,findci=true,
                    bound=true,bound_threshold=0.8,kwds...)
  selection = select_params(params;kwds...)
  Nst = length(unique(params.Δf))
  if length(selection) != Nst
    error("Expected $Nst parameter entires, one for each Δf.",
          "\nInstead found the entires: ",string(selection),
          "\nKeyword selection: ",string(kwds))
  end
  dfstream_ind = @linq df |>
    where(in.(:pindex,Ref(selection))) |>
    by([:pindex,:created],
       st = params[first(:pindex),:].Δf,
       streaming = streamprop(:percepts,:length,bound=bound,
                              threshold=bound_threshold),
       streaming_unbound = streamprop(:percepts,:length,bound=false))

  dfstream = by(dfstream_ind,:st) do dfind
    @assert ismissing(resampling_N) || resampling_N == size(dfind,1)

    if any(!ismissing,dfind[:streaming])
      @with dfind begin
        streaming = collect(skipmissing(:streaming[resampling]))
        if findci
          if length(streaming) > 3 &&
              any(x != streaming[1] for x in streaming)
            cint = dbootconf(streaming)
            DataFrame(mean = mean(streaming),
                      lowerc = cint[1],
                      upperc = cint[2])
          else
            DataFrame(mean = mean(streaming),
                      lowerc = minimum(streaming),
                      upperc = maximum(streaming))
          end
        else
          DataFrame(mean = mean(streaming))
        end
      end
    else
      # if we don't have any (or very little) bistability at all
      # the analysis above will fail: we still want to report the
      # dominant percept (but no estimate of bounds)
      val = @with(dfind,mean(:streaming_unbound[resampling]))
      if findci
        DataFrame(mean = val, lowerc=val, upperc=val)
      else
        DataFrame(mean = val)
      end
    end
  end

  dfstream.experiment = "simulation"

  if keep_simulations
    dfstream, dfstream_ind
  else
    dfstream
  end
end


function stream_dfh(;findci=true)
  dfh = if findci
    @by(CSV.read(joinpath("..","analysis","context","stream_prop.csv")),:st,
        mean = mean(skipmissing(:response)),
        rms = rms(mean(skipmissing(:response)) .- skipmissing(:response)),
    lowerc = dbootconf(collect(skipmissing(:response)))[1],
    upperc = dbootconf(collect(skipmissing(:response)))[2])
  else
    @by(CSV.read(joinpath("..","analysis","context","stream_prop.csv")),:st,
        mean = mean(skipmissing(:response)),
        rms = rms(mean(skipmissing(:response)) .- skipmissing(:response)))
  end

  dfh[:experiment] = "human"
  dfh
end

addtuple(x,y) = NamedTuple{keys(x)}(map(+,x,y))

function stream_rms(df,params,dfh;N=1000,return_parts=false,kwds...)
  dfsm, dfs = stream_dfs(df,params;keep_simulations=true,findci=false,kwds...)
  dfh = @where(dfh,:st .∈ Ref(dfsm.st))

  # shuffle the indices ensures that the measure does not reflect any accidental
  # correlation across stimulus conditions (since we know there is nothing
  # special about the ordering)
  result = mapreduce(addtuple,1:N) do n
    dfs[:index] = 0
    for g in groupby(dfs,:st)
      g.index = shuffle(1:size(g,1))
    end

    means = by(dfs,[:index]) do dfind
      if size(dfind.streaming,1) != size(dfh.mean,1)
        @warn("Found some missing simulations. We got $(dfind), "*
              "but this was expected to have $(size(dfh.mean,1)) rows.")
        DataFrame(var = missing, rms = missing,
                  rms_unbound = missing, var_unbound = missing)
      else
        DataFrame(var = meansqr(dfsm.mean .- dfind.streaming),
                  rms = rms(dfh.mean .- dfind.streaming),
                  rms_unbound = rms(dfh.mean .- dfind.streaming_unbound),
                  var_unbound = meansqr(dfsm.mean .- dfind.streaming_unbound))
      end
    end

    if all(!ismissing,means.rms)
      (rms = mean(skipmissing(means.rms)), var = meansqr(skipmissing(means.var)))
    else
      (rms = mean(means.rms_unbound), var = meansqr(means.var_unbound))
    end
  end

  if return_parts
    (rms = (result.rms/N) / rms(dfh.rms),
     sd = √(result.var/N), #/ rms(dfh.rms),
     mratio = (result.rms/N) / rms(dfsm.mean .- dfh.mean))
  else
    (result.rms/N) / rms(dfh.rms)
  end
end

const N_for_pressnitzer_hupe_2006 = 23
# manually selected range to ensure reasonable bin size given the bin size of
# the original data from P&H 2006
const normalized_hist_range = 0:0.33333:20

function normweights(x)
  total = sum(x)
  total > 0 ? x ./ total : x
end

function length_rms(df,params,human;return_parts = false,kwds...)
  len, unlen = length_sdata(df,params;kwds...)
  hist = fit(Histogram,len,weights(len),normalized_hist_range)

  dens = normweights(hist.weights)
  serr = rms(human.dens .- dens)

  if return_parts
    (rms = serr / human.err, sd = std(log.(unlen)))
  else
    serr / human.err
  end
end

function length_dfh()
  ph = CSV.read(joinpath("..","data","pressnitzer_hupe",
                         "pressnitzer_hupe_inferred.csv"))
  DataFrame(length = ph.length,experiment="human",
            nlength = normlength(ph.length))
end

function length_hdata()
  ph = CSV.read(joinpath("..","data","pressnitzer_hupe",
                         "pressnitzer_hupe_inferred.csv"))
  len = normlength(ph.length)

  hist = fit(Histogram,len,weights(len),normalized_hist_range)
  dens = hist.weights ./ sum(hist.weights)

  nsamples = div(length(len),N_for_pressnitzer_hupe_2006)
  err = mean(dbootinds(len,numobsperresample=nsamples,numresample=10^4)) do inds
    hist = fit(Histogram,len[inds],weights(len[inds]),normalized_hist_range)
    dens_ = normweights(hist.weights)
    rms(dens_ .- dens)
  end

  (dens=dens,err=err)
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

function length_sdata(df,params;resample_N,resampling=Colon(),kwds...)
  selection = select_params(params;Δf=6,kwds...)
  dfs = @where(df,(:pindex .== selection))
  if resampling isa Colon
    normlength(dfs.length), dfs.length
  else
    # TODO:
  end
end

function length_dfs(df,params;kwds...)
  selection = select_params(params;Δf=6,kwds...)
  dfs = @where(df,(:pindex .== selection))
  DataFrame(length = dfs.length,experiment="simulation",
            nlength = normlength(dfs.length));
end

length_df(df,params;kwds...) = vcat(length_dfh(), length_dfs(df,params;kwds...))

function plot_fitmask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                      stop_time=20s,kwds...)
  fit = plot_fit(df,params;kwds...)
  mask = plot_mask(df,params,settings;Δf=Δf,simulation=simulation,
                   start_time=start_time,stop_time=stop_time,kwds...)
  vstack(fit,mask)
end

function handlebound(fn,seconds;bound=true,threshold=0.8)
    if bound && length(seconds) < 3
        return missing
    end

    if !bound || (sum(seconds[2:end-1]) > threshold*sum(seconds))
        fn(1:length(seconds))
    else
        fn(2:length(seconds)-1)
    end
end

function buildup_mean(buildup_df;delta,length)
  buildup = DataFrame(time=range(0,stop=length,step=delta))
  buildup[:value] = 0.0
  runs = groupby(buildup_df,:run)
  for run in runs
    j = 1
    ts = cumsum(run.length)
    for (i,t) in enumerate(buildup.time)
      while j <= Base.length(ts) && t > ts[j]
        j += 1
      end
      j <= Base.length(ts) || break
      buildup.value[i] += (buildup_df.response[j]-1)
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
