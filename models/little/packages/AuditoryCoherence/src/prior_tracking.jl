using Combinatorics

include("tracking_priors.jl")
@with_kw struct PriorTracking <: Tracking
  tc::typeof(1.0s) = 1s
  thresh::Float64 = 1e-2
  source_prior
  freq_prior
  unmodeled_prior::Float64 = 1e-2
  max_sources::Int = 5
end
Tracking(::Val{:prior};params...) = PriorTracking(;params...)

# TODO: if we go back to this we probably want to allow
# for "unmodeled" sources
# ALSO note that i and k were switched when I used this function
# (have to think about what that means below
possible_orders(N,c) =
  (order for selection in combinations(1:N,c)
    for order in permutations(selection))

# the number of ways an improper subset
# of observations A can be mapped to a set of
# source in B when it is possible that
# source b ∈ B is modeled as sum of elements in A
struct Grouping
  groups::Vector{Vector{Int}}
  sources::Vector{Int}
end
possible_groupings(n_sources,n_obs) =
  (Grouping(grouping,mapping)
   for observations in combinations(1:n_obs)
   for grouping in partitions(observations)
   for sources in combinations(1:n_sources,length(grouping))
   for mapping in permutations(sources))

components(x) = vcat(x.groups...)
Base.sum(fn::Function,x::Grouping) = sum(fn,zip(x.groups,x.sources))
iterable(x::Grouping) = zip(x.groups,x.sources)

struct TrackedSources{S,F}
  sources::Vector{S}
  freqs::Vector{F}
  params::PriorTracking
  last_observed::Vector{Bool}
end
function TrackedSources(params::PriorTracking)
  TrackedSources([zero(params.source_prior) for i in 1:params.max_sources],
                 [BinomialCond(:old => Beta(0.0,0.0),:new => Beta(0.0,0.0))
                  for i in 1:params.max_sources],
                 params,
                 [false for i in 1:params.max_sources])
end

function logpdf(track::TrackedSources,C::AxisArray,grouping::Grouping)
  # observed, modeled sources
  logsum = sum(grouping) do kk_i
    (kk,i) = kk_i
    observed = sum(k -> vec(component(C,k)),kk)
    logpdf_thresh(track.params.source_prior+track.sources[i], observed, track.params.thresh) +
    logpdf(track.params.freq_prior + track.freqs[i],
           track.last_observed[i] ? :old : :new,true)
  end

  # unobserved (but modeled) sources
  unobserved = setdiff(1:track.params.max_sources,grouping.sources)
  if !isempty(unobserved)
    logsum += sum(unobserved) do i
      logpdf(track.params.freq_prior + track.freqs[i],
             track.last_observed[i] ? :new : :old,false)
    end
  end

  # unmodeled (but observed) components
  unmodeled = setdiff(1:ncomponents(C),components(grouping))
  if !isempty(unmodeled)
    logsum += sum(unmodeled) do k
      logpdf_thresh(track.params.source_prior,vec(component(C,k)),
                    track.params.thresh)
    end
    logsum += log(track.params.unmodeled_prior)
  end

  logsum
end

function mult!(track::TrackedSources,x)
  for i in 1:track.params.max_sources
    mult!(track.sources[i],x)
    mult!(track.freqs[i],x)
  end
end

function update!(track::TrackedSources,Csource,grouping,w=1.0)
  for i in grouping.sources
    update!(track.sources[i],vec(component(Csource,i)),w)
    update!(track.freqs[i],track.last_observed[i] ?  :old : :new,true,w)
    track.last_observed[i] = true
  end
  for i in setdiff(1:track.params.max_sources,grouping.sources)
    update!(track.freqs[i],track.last_observed[i] ? :new : :old,false,w)
    track.last_observed[i] = false
  end
end


function track(C::Coherence,params::PriorTracking)
  @show params.max_sources
  time = Axis{:time}
  component = Axis{:component}
  # hopefully one day I can remove these assertions
  # but not right now
  @assert axisdim(C,time) == 1
  @assert axisdim(C,component) == 4

  track = TrackedSources(params)

  C_out = similar(C,axes(C)[1:end-1]...,component(1:params.max_sources))
  C_out .= 0
  source_out = copy(C_out)
  sourceS_out = copy(C_out)
  freq_out = zeros(ntimes(C),params.max_sources,4)

  @show decay = 1.0 - 1.0 / max(1.0,params.tc / Δt(C))

  @showprogress "Tracking Sources..." for t in eachindex(times(C))
    # find the MAP grouping of sources
    MAPgrouping = maximumby(grouping -> logpdf(track,C[time(t)],grouping),
                            possible_groupings(params.max_sources,
                                               ncomponents(C)))

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
    mult!(track,decay)
    update!(track,C_out[time(t)],MAPgrouping)

    for i in MAPgrouping.sources
      source_out[t,:,:,i] = track.sources[i].μ
      sourceS_out[t,:,:,i] = std(track.sources[i])

      freq_out[t,i,1] = track.freqs[i].at[:old].α
      freq_out[t,i,2] = track.freqs[i].at[:old].β
      freq_out[t,i,3] = track.freqs[i].at[:new].α
      freq_out[t,i,4] = track.freqs[i].at[:new].β
    end
  end

  # rearrange so largest is first
  order = sortperm(component_means(C_out),rev=true)
  C_out .= C_out[component(order)]

  C_out,source_out,freq_out,sourceS_out
end
