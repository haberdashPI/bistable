export map_components

@with_kw struct MultiPriorTracking <: Tracking
  tcs::Array{typeof(1.0s)} = [100ms, 250ms, 500ms, 1s, 2.5s, 5s]
  source_priors
  freq_prior
  unmodeled_prior::Float64 = 0.0
  max_sources::Int = 4
end
Tracking(::Val{:multi_prior};params...) = MultiPriorTracking(;params...)

function expand_params(params::MultiPriorTracking)
  [PriorTracking(tc,prior,params.freq_prior,
                 params.unmodeled_prior,params.max_sources)
   for tc in params.tcs for prior in params.source_priors]
end

function nitr(C::Coherence,params::MultiPriorTracking)
  ntimes(C) * length(params.tcs) * length(params.source_priors)
end

function track(C::Coherence,params::MultiPriorTracking,progressbar=true,
               progress=track_progress(progressbar,nitr(C,params),"multi-prior"))
  map(expand_params(params)) do p
    track(C,p,true,progress)
  end
end

function map_components(fn,tracks::Array{<:Tuple{Coherence,AbstractArray}};
                        window=500ms,step=250ms)
  Ct1 = tracks[1][1]
  windows = windowing(Ct1,length=window,step=step)
  ts = linspace(times(Ct1)[1],times(Ct1)[end]-step,length(windows))

  result = map(enumerate(windows)) do (i_ixs)
    i,ixs = i_ixs
    best_track = map(tracks) do results
      mean(results[2][ixs])
    end |> indmax

    fn(tracks[best_track][1][Axis{:time}(ixs)])
  end

  AxisArray(result,Axis{:time}(ts))
end

function mask(sp::AuditoryModel.AuditorySpectrogram,
              tracks::Array{<:Tuple{Coherence,AbstractArray}},
              settings;progressbar=false,kwds...)
  scales = settings["scales"]["values"] .* cycoct
  freql,freqh = settings["rates"]["freq_limits"] .* Hz

  cr = cortical(sp,scales=scales,progressbar=progressbar)
  cr = cr[:,:,freql .. freqh]

  mask(cr,tracks;progressbar=progressbar,kwds...)
end

function mask(cr::AuditoryModel.Cortical,
              tracks::Array{<:Tuple{Coherence,AbstractArray}};
              order=1,window=500ms,step=250ms,progressbar=false)

  Ct1 = tracks[1][1]
  @assert axisdim(cr,Axis{:time}) == 1
  @assert axisdim(cr,Axis{:scale}) == 2
  @assert axisdim(cr,Axis{:freq}) == 3
  @assert size(cr)[2:end] == size(Ct1)[2:end-1] "Dimension mismatch"

  windows = windowing(Ct1,length=window,step=step)
  ts = linspace(times(Ct1)[1],times(Ct1)[end]-step,length(windows))

  progress = progressbar ? Progress(length(ts),"Masking: ") : nothing
  mask_helper(cr,tracks,order,windows,ts,progress)
end

function mask_helper(cr,tracks,order,windows,ts,progress)
  y = zeros(AxisArray(cr))
  norm = similar(y,real(eltype(cr)))
  norm .= zero(real(eltype(cr)))

  cohere_windows = collect(windowing(cr,AuditoryModel.Params(tracks[1][1])))

  for (i,ixs) = enumerate(windows)
    best_track = map(tracks) do results
      mean(results[2][ixs])
    end |> indmax
    components = tracks[best_track][1][Axis{:time}(ixs)]
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


