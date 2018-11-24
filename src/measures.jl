using StatsBase: fit, Histogram, weights
using Random

rms(x) = sqrt(mean(x.^2))

function select_params(params;kwds...)
  condition = trues(size(params,1))
  for (var,val) in pairs(kwds)
    condition .&= abs.(params[var] .- val) .<= 0.1
  end
  params[condition,:pindex]
end

function setup_human_data()
  (stream=stream_dfh(),lengths=length_dfh())
end

function model_rms(df,params,dfh;return_parts=false,kwds...)
  ((stream,stream_mean),lengths) =
  (stream_rms(df,params,dfh.stream,mean_v_ind=true;kwds...),
   length_rms(df,params,dfh.lengths;kwds...))

  if return_parts
    (rms = mean((stream,lengths)),
     stream_rms = stream,
     stream_mean_rms_error = stream_mean,
     length_rms = lengths)
  else
    mean((stream,lengths))
  end
end

function stream_df(df,params;findci=true,bound=true,bound_threshold=0.8,kwds...)
  dfs = stream_dfs(df,params;findci=findci,bound=bound,
                   bound_threshold=bound_threshold,kwds...)
  dfh = stream_dfh(findci=findci)
  delete!(dfh,:rms)
  vcat(dfs,dfh)
end

function stream_dfs(df,params;keep_simulations=false,findci=true,
                    bound=true,bound_threshold=0.8,kwds...)
  selection = select_params(params;kwds...)
  if length(selection) != 3
    error("Expected three parameter entires, one for each Δf.",
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
    if any(!ismissing,dfind[:streaming])
      @with dfind begin
        streaming = collect(skipmissing(:streaming))
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
      val = @with(dfind,mean(:streaming_unbound))
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

function stream_rms(df,params,dfh;N=1000,mean_v_ind=false,kwds...)
# use the individual runs of the simulation
# and get the rms of each: then bootstrap it
dfsm, dfs = stream_dfs(df,params;keep_simulations=true,findci=false,kwds...)

# shuffle the indices ensures that the measure does not reflect any accidental
# correlation across stimulus conditions (since we know there is nothing
# special about the ordering)
result = mean(1:N) do n
  dfs[:index] = 0
  for g in groupby(dfs,:st)
    g.index = shuffle(1:size(g,1))
  end

  means = by(dfs,[:index]) do dfind
    DataFrame(mean = mean(dfind.streaming),
              rms = rms(dfh.mean .- dfind.streaming),
              rms_unbound = rms(dfh.mean .- dfind.streaming_unbound))
  end

  if all(!ismissing,means.rms)
    (mean(skipmissing(means.rms))) ./ (rms(dfh.rms))
  else
    mean(means.rms_unbound) ./ rms(dfh.rms)
  end
end

if mean_v_ind
  result, result - mean(dfsm.mean .- dfh.mean) ./ rms(dfh.rms)
else
  result
end
end

const N_for_pressnitzer_hupe_2006 = 23
# manually selected range to ensure reasonable bin size given the bin size of
# the original data from P&H 2006
const normalized_hist_range = 0:0.166666:10

function length_rms(df,params,dfh;kwds...)
  slen = length_dfs(df,params;kwds...).nlength;
  hlen = dfh.nlength;

  shist = fit(Histogram,slen,weights(slen),normalized_hist_range)
  hhist = fit(Histogram,hlen,weights(hlen),normalized_hist_range)

  hdens = hhist.weights ./ sum(hhist.weights)
  sdens = shist.weights ./ sum(shist.weights)
  serr = rms(hdens .- sdens)

  nsamples = div(length(hlen),N_for_pressnitzer_hupe_2006)
  herr = mean(dbootinds(hlen,numobsperresample=nsamples)) do inds
    hhist = fit(Histogram,hlen[inds],weights(hlen[inds]),normalized_hist_range)
    hdens_ = hhist.weights ./ sum(hhist.weights)
    rms(hdens_ .- hdens)
  end
  serr / herr
end

function length_dfh()
  ph = CSV.read(joinpath("..","data","pressnitzer_hupe",
                         "pressnitzer_hupe_inferred.csv"))
  DataFrame(length = ph[:length],experiment="human",
            nlength = exp.(zscore(log.(ph[:length]))));
end

function length_dfs(df,params;kwds...)
  selection = select_params(params;Δf=6,kwds...)
  dfs = @where(df,(:pindex .== selection))
  DataFrame(length = dfs[:length],experiment="simulation",
            nlength = exp.(zscore(log.(dfs[:length]))));
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
