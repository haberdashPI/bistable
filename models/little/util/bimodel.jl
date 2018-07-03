include(joinpath(@__DIR__,"biapply.jl"))
using AuditoryCoherence

function bistable_model(stim_count,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)

  stim = ab(params[:delta_t]/2,params[:delta_t]/2,1,stim_count,
            params[:standard_f],params[:delta_f]) |>
         normpower |> amplify(-10dB)

  spect = audiospect(stim,progressbar=progressbar)
  spectat = apply_bistable(spect,:freqs,params,settings,progressbar=progressbar,
                           intermediate_results=intermediate_results)
  specta = spectat[1]

  scales = cycoct.*settings["scales"]["values"]
  cs = cortical(spect,
                progressbar=progressbar,scales=scales,
                bandonly=settings["config"]["bandonly"])

  csat = apply_bistable(cs,:scales,params,settings,progressbar=progressbar,
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
    results...,tracks,C,csat...,spectat...
  else
    results[1]
  end
end

