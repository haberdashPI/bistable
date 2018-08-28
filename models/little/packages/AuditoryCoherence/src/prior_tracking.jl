using Combinatorics
using DataStructures

struct PriorTracking{Sp,Fp} <: Tracking
  tc::typeof(1.0s)
  source_prior::Sp
  freq_prior::Fp
  max_sources::Int
end
function Tracking(C,::Val{:prior};tc=1s,source_prior=nothing,
                  freq_prior=nothing, max_sources=4)
  PriorTracking(tc,source_prior,freq_prior,max_sources)
end

include("tracking_priors.jl")

# the number of ways an improper subset
# of observations A can be mapped to a set of
# source in B when it is possible that
# source b ∈ B is modeled as sum of elements in A
possible_groupings(n_sources,n_obs) =
  (Grouping(grouping,mapping)
   for grouping in partitions(1:n_obs)
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

  track = TrackedSources(prod(size(C,2,3)),params)
  C_ = permutedims(C,[2,3,4,1])
  timeax = 4

  C_out = similar(C_)
  C_out .= 0
  source_out = copy(C_out)
  sourceS_out = copy(C_out)

  window = ceil(Int,params.tc / Δt(C))
  lp_out = AxisArray(fill(0.0,ntimes(C)),axes(C,1))
  buf = Array{eltype(C_)}(size(C_,1,2)...)
  function sumcomponents(C_t,kk)
    y = sum!(buf,view(C_t,:,:,kk,:))
    # buf ./= size(C_t,timeax)
    vec(buf)
  end

  veccomponent(C_t,k) = vec(view(C_t,:,:,k))

  groupings = CircularDeque{Grouping}(window+1)
  for t in eachindex(times(C))
    # find the MAP grouping of sources
    MAPgrouping, lp_out[t] = # OLD PROFILE COUNT: 18922
      maximumby(g -> logpdf(track,view(C_,:,:,:,t:t),sumcomponents,g),
                possible_groupings(params.max_sources,ncomponents(C)))

    # arrange the sources
    for (kk,i) in iterable(MAPgrouping)
      for k in kk
        C_out[:,:,i,t] .+= C_[:,:,k,t]  # PROFILE COUNT: 2040
      end
    end

    # update the source models
    update!(track,view(C_out,:,:,:,t),veccomponent,MAPgrouping)
    # mult!(track,decay)
    push!(groupings,MAPgrouping)
    if t > window
      downdate!(track,view(C_out,:,:,:,t-window),veccomponent,shift!(groupings))
    end

    next!(progress)
  end


  C_out_perm = similar(C)
  permutedims!(C_out_perm,C_out,[4,1,2,3])

  (C_out_perm,lp_out)
end

function sort_components(x::Coherence)
  order = sortperm(component_means(x),rev=true) # OLD PROFILE COUNT: 2313
  x .= x[component(order)]
  x
end
