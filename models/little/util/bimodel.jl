include(joinpath(@__DIR__,"biscales.jl"))
using AuditoryCoherence

function count_streams(tracks;window=500ms,step=250ms,threshold=2,min_length=1s,
                       progressbar=false,intermediate_results=false)
  Ct1 = tracks[1][1]
  @assert axisdim(Ct1,Axis{:time}) == 1
  windows = windowing(Ct1,length=window,step=step)
  ts = linspace(times(Ct1)[1],times(Ct1)[end]-step,length(windows))
  ratios = fill(0.0,length(windows))

  for (i,ixs) in enumerate(windows)
    best_track = map(tracks) do results
      mean(results[2][ixs])
    end |> indmax
    track_window = tracks[best_track][1][Axis{:time}(ixs)]

    strengths = sort(component_means(track_window),rev=true)
    ratios[i] = strengths[1] / sum(strengths[2:end])

  end

  counts = percept_lengths(AxisArray(ratios .> threshold,Axis{:time}(ts)),
                           min_length)
  if intermediate_results
    counts,AxisArray(ratios,Axis{:time}(ts))
  else
    (counts,)
  end
end

function bistable_model(stim_count,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)

  scales = cycoct.*settings["scales"]["values"]
  stim = ab(params[:delta_t]/2,params[:delta_t]/2,1,stim_count,
            params[:standard_f],params[:delta_f]) |>
         normpower |> amplify(-10dB)
  cs = cortical(audiospect(stim,progressbar=progressbar),
                progressbar=progressbar,scales=scales,
                bandonly=settings["config"]["bandonly"])

  csat = bistable_scales(cs,params,settings,progressbar=progressbar,
                        intermediate_results=intermediate_results)
  csa = csat[1]
  rates = Float64.(settings["rates"]["values"]).*Hz
  start,stop = settings["rates"]["freq_limits"]
  crs = cortical(csa[:,:,start*Hz .. stop*Hz], rates=[-rates;rates],
                 bandonly=settings["config"]["bandonly"],
                 progressbar=progressbar)

  C = cohere(
    crs,
    ncomponents=settings["nmf"]["K"],
    window=settings["nmf"]["window_ms"]*ms,
    method=:nmf,
    skipframes=settings["nmf"]["skipframes"],
    delta=settings["nmf"]["delta_ms"]*ms,
    maxiter=settings["nmf"]["maxiter"],
    tol=settings["nmf"]["tol"],
    progressbar=progressbar
  )

  tracks = track(
    C,
    method=:multi_prior,
    max_sources = settings["track"]["max_sources"],
    tcs = settings["track"]["time_constants_s"].*s,

    source_priors = [isonorm(sd,N,size(C,2,3))
                     for sd in settings["track"]["source_prior"]["sds"]
                     for N in settings["track"]["source_prior"]["Ns"]],

    freq_prior = freqprior(settings["track"]["freq_prior"]["bias"],
                           settings["track"]["freq_prior"]["N"]),
    progressbar=progressbar
  )

  results = count_streams(
    tracks,
    window=settings["percept_lengths"]["window_ms"].*ms,
    step=settings["percept_lengths"]["delta_ms"].*ms,
    threshold=settings["percept_lengths"]["threshold"],
    min_length=settings["percept_lengths"]["min_length_ms"].*ms,
    intermediate_results=intermediate_results,
  )

  if intermediate_results
    results...,tracks,C,csat...
  else
    results[1]
  end
end

