export track
using Combinatorics

abstract type Tracking end
@with_kw struct SimpleTracking <: Tracking
  tc::typeof(1.0s) = 1s
end
Tracking(::Val{:simple};params...) = SimpleTracking(;params...)

function track(C::Coherence;method=:simple,params...)
  method = Tracking(Val{method}();params...)
  track(C,method)
end

function maximumby(by,xs)
  state = start(xs)
  # NOTE: in general shoud check for empty iterator
  # but I know there will be at least one value in usage below
  # and this avoids poor performance in julia 0.6 due to type stability
  # (should be fixed in 0.7)
  # if done(xs,state)
  # return nothing
  # else
  result, state = next(xs,state)
  maxval = by(result)
  while !done(xs,state)
    x, state = next(xs,state)
    val = by(x)
    if isless(maxval,val)
      maxval = val
      result = x
    end
  end
  result
  # end
end

function bestordering(x,y)
  @assert length(x) == length(y)
  @assert length(x) <= 6 "No more than 6 components supported"*
    " (due to combinatorial explosion of possible solutions.)"
  # minimize cosine difference
  diff(x,y) = norm(vec(x) .- vec(y),1)

  maximumby(permutations(1:length(y))) do order
    -sum(diff(xi,yi) for (xi,yi) in zip(x,y[order]))
  end
end

function track(C::Coherence,params::SimpleTracking)
  C = copy(C)
  time = Axis{:time}
  component = Axis{:component}
  K = ncomponents(C)
  λ = Δt(C)/params.tc
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

