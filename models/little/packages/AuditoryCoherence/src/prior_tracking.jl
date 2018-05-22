using Combinatorics

include("tracking_priors.jl")
# TODO: make the prior type more generic (can any of the methods above be
# generic)
@with_kw struct PriorTracking <: Tracking
  tc::typeof(1.0s) = 1s
  thresh::Float64 = 1e-2
  prior
end
Tracking(::Val{:prior};params...) = PriorTracking(;params...)

possible_orders(N,c) =
  (order for selection in combinations(1:N,c)
    for order in permutations(selection))

const max_sources = 5
function track(C::Coherence,params::PriorTracking)
  # TODO: how to deal with the sparse, very near zero
  # values, leading to poorly defined covariance
  # (some sort of sparisty assumption)

  time = Axis{:time}
  component = Axis{:component}
  sources = [zero(params.prior) for i in 1:max_sources]
  freqs = fill(0.0,max_sources)
  count = 0.0
  freq_prior = 1,2

  C_out = similar(C,axes(C)[1:end-1]...,component(1:max_sources))

  @showprogress "Tracking Sources..." for t in eachindex(times(C))
    # find the MAP ordering of sources
    bestorder = maximumby(possible_orders(max_sources,ncomponents(C))) do order
      logsum = sum(enumerate(order)) do i_k
        (i,k) = i_k
        logpdf_thresh(params.prior + sources[k],vec(C[time(t),component(i)]),
                      params.thresh) +
          log((freqs[k] + freq_prior[1])/(count + freq_prior[2]))
      end
      logsum += sum(setdiff(1:max_sources,order)) do i
        log(1-((freqs[i] + freq_prior[1])/(count + freq_prior[2])))
      end

      logsum
    end

    # rearrange the sources
    C_out[time(t),component(bestorder)] = C[time(t)]
    C_out[time(t),component(setdiff(1:max_sources,bestorder))] = 0

    tc = 1 / max(0.5,params.tc / Δt(C))
    # update the source models

    for k in 1:max_sources
      mult!(sources[k],1-tc)
      freqs[k] *= (1-tc)
    end
    count *= (1-tc)

    for (i,k) in enumerate(bestorder)
      update!(sources[k],vec(C[time(t),component(i)]))
      freqs[k] += tc
    end
    count += tc
  end

  # rearrange so largest is first
  ordering = sortperm(component_means(C_out),rev=true)
  C_out .= C_out[component(ordering)]
end
