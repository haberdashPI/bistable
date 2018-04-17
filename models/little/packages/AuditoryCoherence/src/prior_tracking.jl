
using Combinatorics

mutable struct MultiNormalStats{T}
  x::Vector{T}
  x2::Matrix{T}
  n::Float64
  x2_offset::Float64
end
function MultiNormalStats(data::AbstractMatrix,n=size(data,1),
                          x2_offset=size(data,2))
  x = squeeze(sum(data,1),1) .* (n/size(data,1))
  x2 = (data'data) .* (n/size(data,1))
  MultiNormalStats{eltype(data)}(x,x2,n,x2_offset)
end
MultiNormalStats(::Type{T},n) where T =
  MultiNormalStats{T}(zeros(T,n),zeros(T,n,n),0,0)

function update!(stats::MultiNormalStats{T},x::AbstractVector{T}) where T
  stats.x .+= x
  stats.x2 .+= x.*x'
  stats.n += 1

  stats
end
update(stats::MultiNormalStats{T},x::AbstractVector{T}) where T =
  update!(copy(stats),x)

function mult!(stats::MultiNormalStats,c)
  stats.x *= c
  stats.x2 .*= c
  stats.n *= c

  stats
end

Base.:(+)(a::MultiNormalStats{T},b::MultiNormalStats{T}) where T =
  MultiNormalStats(a.x+b.x,a.x2+b.x2,a.n+b.n,a.x2_offset+b.x2_offset)

function logpdf(stats::MultiNormalStats,x::AbstractVector)
  d = length(x)
  Λ = (stats.x2 .- (stats.x.*stats.x')) ./ stats.n
  μ = stats.x ./ stats.n
  n = stats.n
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n + stats.x2_offset - d + 1))

  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x)
end
pdf(stats::MultiNormalStats,x::AbstractVector) = exp(logpdf(stats,x))

function logpdf_mvt(v,μ,Σ,x)
  d = length(μ)
  C = lgamma((v+1)/2) - (lgamma(v/2)+log(v*π)^(d/2)) - 0.5logabsdet(Σ)[1]
  diff = abs.(μ.-x)

  C*log(1+1/v*(diff'/Σ)*diff)*-(v+d)/2
end

@with_kw struct PriorTracking <: Tracking
  tc::typeof(1.0s) = 1s
  prior::MultiNormalStats{Float64}
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
  sources = [MultiNormalStats(Float64,length(params.prior.x))
             for i in 1:max_sources]
  freqs = fill(0.0,max_sources)
  count = 0.0
  freq_prior = 1,2

  C_out = similar(C,axes(C)[1:end-1]...,component(1:max_sources))

  @showprogress "Tracking Sources..." for t in eachindex(times(C))
    # find the MAP ordering of sources
    bestorder = maximumby(possible_orders(max_sources,ncomponents(C))) do order
      logsum = sum(enumerate(order)) do i_k
        (i,k) = i_k
        logpdf(params.prior + sources[k],vec(C[time(t),component(i)])) +
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
    @show tc
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
