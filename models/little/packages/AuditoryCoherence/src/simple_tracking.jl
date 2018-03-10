export track

function orderings(N)
  if N > 5
    error("Orderings larger than 5 not supported")
  end
  insert!(array,x,i) = [j == i ? x : array[j < i ? j : j-1]
                        for j in 1:length(array)+1]
  if N == 0; [Int[]]
  else [insert!(copy(o),N,i) for i in 1:N for o in orderings(N-1)] end
end

function bestordering(x,y)
  @assert length(x) == length(y)
  score(x,y,ordering) = norm(x .- y[ordering],1)
  os = orderings(length(y))
  i = indmin(score(x,y,o) for o in os)

  os[i]
end

function track(C::Coherence;tc=1s)
  time = Axis{:time}
  component = Axis{:component}
  K = ncomponents(C)
  λ = Δt(C)/tc
  Ĉ = C[time(1)]

  @showprogress "Tracking Sources..." for t in eachindex(times(C))
    ordering = bestordering([Ĉ[component(k)] for k in 1:K],
                            [C[time(t),component(k)] for k in 1:K])
    C[time(t)] = C[time(t),component(ordering)]
    Ĉ .= (1-λ).*Ĉ .+ (λ).*C[time(t)]
  end

  # rearrange so largest is first

  ordering = sortperm(component_means(C),rev=true)
  C .= C[component(ordering)]
end

# next simplest
# steps:

# at each step take an MAP approach
# 1. for each scale:
#   a. eliminate components below a threshold
#   b. group remaing observed components with their closest source
#   c. any source not included is marked as inactive
#   d. use the model prior to compute a score
# 2. select the model with the highest score

# NOTE: not quite gonna work, since NMF doens't work well when there's
# noise.

# function track_sources(C::NMFSeries,params::SimpleTracker)
# end
