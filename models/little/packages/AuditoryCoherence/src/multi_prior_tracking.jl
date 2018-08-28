export map_components

@with_kw struct MultiPriorTracking <: Tracking
  time_constants::Array{typeof(1.0s)}
  source_priors::AxisArray
  freq_prior
  max_sources::Int = 4
end
function Tracking(C,::Val{:multi_prior};time_constants_s=[4],
                  time_constants=time_constants_s*s,
                  source_prior_sds=nothing,source_prior_N=nothing,
                  source_priors=nothing,
                  freq_prior=nothing,
                  freq_prior_N = 2, freq_prior_bias = 0,
                  params...)
  if source_priors == nothing
    @assert(source_prior_sds != nothing,
            "Missing keyword argument `source_prior_sds`.")
    @assert(source_prior_N != nothing,
            "Missing keyword argument `source_prior_N`.")
    source_priors = AxisArray([isonorm(sd,source_prior_N,size(C,2,3))
                               for sd in source_prior_sds],
                              Axis{:prior}(source_prior_sds))
  end

  if freq_prior == nothing
    freq_prior = freqprior(freq_prior_bias,freq_prior_N)
  end
  MultiPriorTracking(;source_priors=source_priors,freq_prior=freq_prior,
                     time_constants=time_constants,params...)
end

function expand_params(params::MultiPriorTracking)
  AxisArray([PriorTracking(tc,prior,params.freq_prior,params.max_sources)
             for tc in params.time_constants for prior in params.source_priors],
            Axis{:params}([(tc,prior) for tc in params.time_constants
                           for prior in axisvalues(params.source_priors)[1]]))
end

function nitr(C::Coherence,params::MultiPriorTracking)
  ntimes(C) * length(params.time_constants) * length(params.source_priors)
end

function track(C::Coherence,params::MultiPriorTracking,progressbar=true,
               progress=track_progress(progressbar,nitr(C,params),"multi-prior"))
  all_params = expand_params(params)
  S = Array{typeof(C)}(size(all_params,1))
  lp = Array{Array{Float64}}(size(all_params,1))
  #=@Threads.threads=# for (i,p) in collect(enumerate(all_params))
    S[i], lp[i] = track(C,p,true,nothing)
  end

  (AxisArray(S, axes(all_params,1)),
   AxisArray(hcat(lp...), axes(C,1), axes(all_params,1)))

end

function map_components(fn,tracks::AxisArray{<:Coherence},
                        tracks_lp::AxisArray{<:Float64};
                        window=500ms,step=250ms)
  # @show size(tracks)
  # @show size(tracks_lp)
  windows = windowing(tracks[1],length=window,step=step)

  result = map(enumerate(windows)) do (i_ixs)
    i,ixs = i_ixs
    best_track = indmax(mean(Array(tracks_lp[ixs,:]),1))
    fn(tracks[best_track][Axis{:time}(ixs)])
  end

  AxisArray(result,axes(windows,Axis{:time}))
end

function mask(sp::AuditoryModel.AuditorySpectrogram,
              tracks::AxisArray{<:Coherence},
              tracks_lp::AxisArray{<:Float64},
              settings;progressbar=false,kwds...)
  scales = settings["scales"]["values"] .* cycoct
  freql,freqh = settings["rates"]["freq_limits"] .* Hz

  cr = cortical(sp,scales=scales,progressbar=progressbar)
  cr = cr[:,:,freql .. freqh]

  mask(cr,tracks,tracks_lp;progressbar=progressbar,kwds...)
end

function mask(cr::AuditoryModel.Cortical,
              tracks::AxisArray{<:Coherence},
              tracks_lp::AxisArray{<:Float64},
              order=1;window=500ms,step=250ms,progressbar=false)

  @assert axisdim(cr,Axis{:time}) == 1
  @assert axisdim(cr,Axis{:scale}) == 2
  @assert axisdim(cr,Axis{:freq}) == 3
  @assert size(cr)[2:end] == size(tracks[1])[2:end-1] "Dimension mismatch"

  windows = windowing(tracks[1],length=window,step=step)

  progress = progressbar ? Progress(length(windows),"Masking: ") : nothing
  mask_helper(cr,tracks,tracks_lp,order,windows,progress)
end

function mask_helper(cr,tracks,tracks_lp,order,windows,progress)
  y = zeros(AxisArray(cr))
  norm = similar(y,real(eltype(cr)))
  norm .= zero(real(eltype(cr)))

  cohere_windows = collect(windowing(cr,AuditoryModel.Params(tracks[1])))

  for (i,ixs) = enumerate(windows)
    best_track = indmax(mean(Array(tracks_lp[ixs,:]),1))
    components = tracks[best_track][Axis{:time}(ixs)]
    sorting = sortperm(component_means(components),rev=true)
    component = components[Axis{:component}(sorting[order])]

    for (ti,t) in enumerate(ixs)
      resh = reshape(component[ti,:,:,:],1,
                     size(component,2:ndims(components)-1...)...)
      y[Axis{:time}(cohere_windows[t])] .+= resh
      norm[Axis{:time}(cohere_windows[t])] += 1
    end

    next!(progress)
  end
  y ./= max.(1,norm)
  y ./= maximum(abs,y)
  y .= (abs.(cr) .* y) .* exp.(angle.(cr).*im)

  cortical(y,AuditoryModel.Params(cr))
end


