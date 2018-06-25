
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
  map(expand_params(params)) do params
    track(C,params,true,progress)
  end
end
