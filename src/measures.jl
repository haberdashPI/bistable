using StatsBase: fit, Histogram, weights
using Random

rms(x) = sqrt(mean(x.^2))


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

