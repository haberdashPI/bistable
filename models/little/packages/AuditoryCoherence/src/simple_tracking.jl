export track, topN
using Combinatorics
using DataStructures

abstract type Tracking end
@with_kw struct SimpleTracking <: Tracking
  tc::typeof(1.0s) = 1s
end
Tracking(::Val{:simple};params...) = SimpleTracking(;params...)

function track(C::Coherence;method=:simple,progressbar=true,params...)
  method = Tracking(Val{method}();params...)
  track(C,method,progressbar)
end

function topN(by,N,xs)
  state = start(xs)
  if done(xs,state)
    error("Empty iterable.")
  end
  x, state = next(xs,state)
  by_x = by(x)
  queue = PriorityQueue{typeof(by_x),typeof(x)}()
  enqueue!(queue,x,by_x)

  while !done(xs,state)
    x, state = next(xs,state)
    by_x = by(x)
    enqueue!(queue,x,by_x)
    if length(queue) > N
      dequeue!(queue)
    end
  end

  result = Array{eltype(xs)}(length(queue))
  i = 0
  while !isempty(queue)
    result[i += 1] = dequeue!(queue)
  end
  result
end

function maximumby(by,xs)
  state = start(xs)
  if done(xs,state)
    error("Empty iterable.")
  end

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

function track_progress(progressbar,n,name)
  if progressbar
    Progress(n,desc="Source Tracking ($name): ")
  end
end

function track(C::Coherence,params::SimpleTracking,progressbar=true,
               progress = track_progress(progressbar,ntimes(C),"simple"))
  C = copy(C)
  time = Axis{:time}
  component = Axis{:component}
  K = ncomponents(C)
  λ = Δt(C)/params.tc
  Ĉ = C[time(1)]

  for t in eachindex(times(C))
    ordering = bestordering([Ĉ[component(k)] for k in 1:K],
                            [C[time(t),component(k)] for k in 1:K])
    C[time(t)] = C[time(t),component(ordering)]
    Ĉ .= (1-λ).*Ĉ .+ (λ).*C[time(t)]

    next!(progress)
  end

  # rearrange so largest is first
  ordering = sortperm(component_means(C),rev=true)
  C .= C[component(ordering)]
end

