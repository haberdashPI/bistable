export bistable_model

function bistable_model(stim_count::Int,params,settings;kwds...)
  # stimulus generation
  stim = ab(params[:Δt]/2,params[:Δt]/2,1,stim_count,params[:f],params[:Δf]) |>
         normpower |> amplify(-10dB)
  bistable_model(stim,params,settings;kwds...)
end

function tokwds(x::Dict{String,<:Any})
  (Symbol(k) => v for (k,v) in x)
end

function bistable_model(stim::AbstractArray,params,settings;interactive=false,
                        progressbar=interactive,
                        intermediate_results=interactive)
  @assert params[:condition] ∈ [:freqs,:scales,:track,:none,:scales_track]

  # auditory spectrogram
  # TODO: configure framelength
  spect = audiospect(stim,progressbar=progressbar;
                     tokwds(settings["freqs"]["analyze"])...)
  spectat = apply_bistable(spect,:freqs,params,progressbar=progressbar,
                           intermediate_results=intermediate_results;
                           tokwds(settings["freqs"]["bistable"])...)
  specta = spectat[1]

  # cortical scales
  cs = cortical(specta, progressbar=progressbar;
                tokwds(settings["scales"]["analyze"])...)
  csclean = cortical(spect, progressbar=progressbar;
                     tokwds(settings["scales"]["analyze"])...)
  csat = apply_bistable(cs,:scales,params,progressbar=progressbar,
                        intermediate_results=intermediate_results;
                        tokwds(settings["scales"]["bistable"])...)

  # cortical rates
  csa = csat[1]
  crs = cortical(csa, progressbar=progressbar; tokwds(settings["rates"])...)

  # temporal coherence (simultaneous grouping)
  C = cohere(crs, method=:nmf, progressbar=progressbar;
             tokwds(settings["nmf"])...)

  # source tracking (sequential grouping)
  tracks,track_lp = track(C, method=:multi_prior, progressbar=progressbar;
                          tokwds(settings["track"]["analyze"])...)

  if params[:condition] == :track || params[:condition] == :scales_track
    after_buildup = settings["track"]["buildup_time_s"]*s .. last(times(C))
    track_lp .-= minimum(track_lp[after_buildup])
    track_lp ./= maximum(track_lp[after_buildup])

    track_lp_at = apply_bistable(track_lp,:track,params,
                              progressbar=progressbar,
                              intermediate_results=intermediate_results;
                              tokwds(settings["track"]["bistable"])...)
    track_lp = track_lp_at[1]
  else
    track_lp_at = (track_lp,)
  end

  # decision making
  ratio = component_ratio(
    tracks,track_lp,
    window=settings["percept_lengths"]["window_ms"].*ms,
    step=settings["percept_lengths"]["delta_ms"].*ms,
    intermediate_results=intermediate_results,
  )

  startHz, stopHz = settings["rates"]["freq_limits_Hz"].*Hz
  bratio = component_bandwidth_ratio(
    csclean[:,:,startHz .. stopHz],
    tracks,track_lp,
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
    counts,ratio,bratio,tracks,track_lp_at...,C,csat...,spectat...
  else
    counts,ratio,bratio
  end
end

