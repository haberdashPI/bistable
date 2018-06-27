struct Grouping
  groups::Vector{Vector{Int}}
  sources::Vector{Int}
end
components(x) = vcat(x.groups...)
Base.sum(fn::Function,x::Grouping) = sum(fn,zip(x.groups,x.sources))
iterable(x::Grouping) = zip(x.groups,x.sources)

struct TrackedSources{S,F}
  sources::Vector{S}
  freqs::Vector{F}
  params::PriorTracking
end
function TrackedSources(C::Coherence,params::PriorTracking)
  TrackedSources([zero(params.source_prior,C) for i in 1:params.max_sources],
                 [Beta(0.0,0.0) for i in 1:params.max_sources],params)
end

function logpdf(track::TrackedSources,C::AxisArray,grouping::Grouping)
  # observed, modeled sources
  logsum = sum(grouping) do kk_i
    (kk,i) = kk_i
    observed = sum(k -> vec(component(C,k)),kk)
    logpdf(track.params.source_prior+track.sources[i], observed) +
      logpdf(track.params.freq_prior + track.freqs[i],true)
  end

  # unobserved (but modeled) sources
  unobserved = setdiff(1:track.params.max_sources,grouping.sources)
  if !isempty(unobserved)
    logsum += sum(unobserved) do i
      logpdf(track.params.freq_prior + track.freqs[i],false)
    end
  end

  # unmodeled (but observed) components
  unmodeled = setdiff(1:ncomponents(C),components(grouping))
  if !isempty(unmodeled)
    logsum += sum(unmodeled) do k
      logpdf(track.params.source_prior,vec(component(C,k)))
    end
    logsum += log(track.params.unmodeled_prior)
  end

  logsum
end

function mult!(track::TrackedSources,x)
  for i in 1:track.params.max_sources
    mult!(track.sources[i],x)
    track.freqs[i] *= x
  end
end

function update!(track::TrackedSources,Csource,grouping,w=1.0)
  for i in grouping.sources
    update!(track.sources[i],vec(component(Csource,i)),w)
    track.freqs[i] = update(track.freqs[i],true,w)
  end
  for i in setdiff(1:track.params.max_sources,grouping.sources)
    track.freqs[i] = update(track.freqs[i],false,w)
  end
end

function downdate!(track::TrackedSources,Csource,grouping,w=1.0)
  for i in grouping.sources
    downdate!(track.sources[i],vec(component(Csource,i)),w)
    track.freqs[i] = downdate(track.freqs[i],true,w)
  end
  for i in setdiff(1:track.params.max_sources,grouping.sources)
    track.freqs[i] = downdate(track.freqs[i],false,w)
  end
end
