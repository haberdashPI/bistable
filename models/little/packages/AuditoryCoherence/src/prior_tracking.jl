using Combinatorics
using DataStructures

@with_kw struct PriorTracking <: Tracking
  tc::typeof(1.0s) = 1s
  source_prior
  freq_prior
  unmodeled_prior::Float64 = 0.0
  max_sources::Int = 4
end
Tracking(::Val{:prior};params...) = PriorTracking(;params...)

include("tracking_priors.jl")

# the number of ways an improper subset
# of observations A can be mapped to a set of
# source in B when it is possible that
# source b ∈ B is modeled as sum of elements in A
possible_groupings(n_sources,n_obs) =
  (Grouping(grouping,mapping)
   for observations in combinations(1:n_obs)
   for grouping in partitions(observations)
   for sources in combinations(1:n_sources,length(grouping))
   for mapping in permutations(sources))

function track(C::Coherence,params::PriorTracking,progressbar=true,
               progress = track_progress(progressbar,ntimes(C),"prior"))
  time = Axis{:time}
  component = Axis{:component}
  # hopefully one day I can remove these assertions
  # but not right now
  @assert axisdim(C,time) == 1
  @assert axisdim(C,component) == 4

  track = TrackedSources(C,params)

  C_out = similar(C,axes(C)[1:end-1]...,component(1:params.max_sources))
  C_out .= 0
  source_out = copy(C_out)
  sourceS_out = copy(C_out)

  window = ceil(Int,params.tc / Δt(C))
  lp_out = AxisArray(fill(0.0,ntimes(C)),axes(C,1))

  groupings = CircularDeque{Grouping}(window+1)
  for t in eachindex(times(C))
    # find the MAP grouping of sources
    MAPgrouping = maximumby(grouping -> logpdf(track,C[time(t)],grouping),
                            possible_groupings(params.max_sources,
                                               ncomponents(C)))
    lp_out[time(t)] = logpdf(track,C[time(t)],MAPgrouping)

    # arrange the sources
    for (kk,i) in iterable(MAPgrouping)
      for k in kk
        #  C_out[time(t),component(i)] .+= C[time(t),component(k)]
        # I hope in the future I can use the above but not now
        # TODO: submit an issue with AxisArrays
        C_out[t,:,:,i] .+= C[t,:,:,k]
      end
    end

    # update the source models
    update!(track,C_out[time(t)],MAPgrouping)
    push!(groupings,MAPgrouping)
    if t > window
      downdate!(track,C_out[time(t-window)],shift!(groupings))
    end

    for i in 1:params.max_sources
      source_out[t,:,:,i] = (track.params.source_prior + track.sources[i]).μ
      sourceS_out[t,:,:,i] = std(track.params.source_prior + track.sources[i])
    end

    next!(progress)
  end

  # rearrange so largest is first
  order = sortperm(component_means(C_out),rev=true)
  C_out .= C_out[component(order)]

  C_out,lp_out,source_out,sourceS_out,track
end
