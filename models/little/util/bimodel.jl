include(joinpath(@__DIR__,"biapply.jl"))
using AuditoryCoherence

function bistable_model(stim_count,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)
  @assert params[:condition] ∈ [:freqs,:scales,:track,:none,:scales_track]

  # stimulus generation
  stim = ab(params[:Δt]/2,params[:Δt]/2,1,stim_count,params[:f],params[:Δf]) |>
         normpower |> amplify(-10dB)

  # auditory spectrogram
  spect = audiospect(stim,progressbar=progressbar)
  spectat = apply_bistable(spect,:freqs,params,settings,progressbar=progressbar,
                           intermediate_results=intermediate_results)
  specta = spectat[1]

  # cortical scales
  scales = cycoct.*settings["scales"]["values"]
  cs = cortical(specta,
                progressbar=progressbar,scales=scales,
                bandonly=settings["config"]["bandonly"])
  csclean = cortical(spect,
                     progressbar=progressbar,scales=scales,
                     bandonly=settings["config"]["bandonly"])
  csat = apply_bistable(cs,:scales,params,settings,progressbar=progressbar,
                        intermediate_results=intermediate_results)

  # cortical rates
  csa = csat[1]
  rates = Float64.(settings["rates"]["values"]).*Hz
  startHz,stopHz = settings["rates"]["freq_limits"].*Hz
  crs = cortical(csa[:,:,startHz .. stopHz], rates=[-rates;rates],
                 bandonly=settings["config"]["bandonly"],
                 progressbar=progressbar)

  # temporal coherence (simultaneous grouping)
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

  # source tracking (sequential grouping)
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

  if params[:condition] == :track || params[:condition] == :scales_track
    track_lp = AxisArray(hcat((x[2] for x in tracks)...),
                         Axis{:time}(times(C)),
                         Axis{:prior}(settings["track"]["source_prior"]["sds"]))
    after_buildup = settings["track"]["buildup_time_s"]*s .. last(times(C))
    track_lp .-= minimum(track_lp[after_buildup])
    track_lp ./= maximum(track_lp[after_buildup])

    track_lp_at = apply_bistable(track_lp,:track,params,settings,
                                 progressbar=progressbar,
                                 intermediate_results=intermediate_results)
    tracksa = collect(zip((x[1] for x in tracks),
                          (view(track_lp_at[1],:,i)
                           for i in 1:size(track_lp_at[1],2))))
  else
    tracksa = tracks
    track_lp_at = (nothing,)
  end

  # decision making
  ratio = component_ratio(
    tracksa,
    window=settings["percept_lengths"]["window_ms"].*ms,
    step=settings["percept_lengths"]["delta_ms"].*ms,
    intermediate_results=intermediate_results,
  )

  bratio = component_bandwidth_ratio(
    csclean[:,:,startHz .. stopHz],
    tracksa,
    window=settings["percept_lengths"]["window_ms"].*ms,
    step=settings["percept_lengths"]["delta_ms"].*ms,
    threshold=settings["percept_lengths"]["bandwidth_threshold"]
  )

  threshold = settings["percept_lengths"]["threshold"]
  # the counts are not perfect at this point but they are used to diagnose
  # serious problems in a simulation, subsequent analysis will examine `ratio`
  # and `bratio` across various thresholds
  counts = percept_lengths(AxisArray(ratio .< threshold,
                                     axes(ratio,Axis{:time})),
                           settings["percept_lengths"]["min_length_ms"].*ms)

  if intermediate_results
    counts,ratio,bratio,tracksa,track_lp_at[2:end]...,C,csat...,spectat...
  else
    counts,ratio,bratio
  end
end

