export isonorm

mutable struct IsoMultiNormalStats{T} <: Stats{T}
  μ::Vector{T}
  S::Vector{T}
  n1::Float64
  n2::Float64
  α::Float64
end
Base.std(x::IsoMultiNormalStats) = x.n2 > 0 ? x.S ./ x.n2 : Inf
Base.zero(x::IsoMultiNormalStats{T},C::Coherence) where T =
  IsoMultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0,0.0)

function isonorm(prior::Coherence,n,α=1,thresh=1e-3)
  data = reshape(mean(prior,4),size(prior,1),:)
  μ = squeeze(mean(data,1),1)
  S = squeeze(sum(data.^2,1),1) .* n./size(data,1) .+ thresh
  IsoMultiNormalStats{eltype(data)}(μ,S,n,n,α)
end

function update!(stats::IsoMultiNormalStats{T},x::AbstractVector{T},w=1.0) where T
  δ = x .- stats.μ
  stats.n1 += w
  stats.μ .+= δ * (w/stats.n1)

  δ2 = x .- stats.μ
  stats.n2 += w^2
  stats.S .+= δ.*δ2 * (w^2/stats.n2)

  stats
end

function mult!(stats::IsoMultiNormalStats,c)
  stats.n1 *= c
  stats.n2 *= c^2

  stats
end

function Base.:(+)(a::IsoMultiNormalStats{T},b::IsoMultiNormalStats{T}) where T
  n1 = a.n1 + b.n1
  n2 = a.n2 + b.n2
  IsoMultiNormalStats(a.μ.*(a.n1./n1) .+ b.μ.*(b.n1./n1), a.S.+b.S,
                      n1,n2,a.α+b.α)
end

function logpdf(stats::IsoMultiNormalStats,x::AbstractVector)
  μ = stats.μ
  α = stats.n1/2 + stats.α
  β = 0.5stats.S
  σ = β./α .* ((stats.n1 + 1) ./ stats.n1)

  sum(logpdf(TDist(2α),(x .- μ)./σ))
end

struct ConstIsoPrior{T} <: Stats{T}
  S::T
  N::Float64
  α::Float64
end
function isonorm(prior::Number,N,α=1,thresh=1e-3)
  ConstIsoPrior{typeof(prior)}(prior.+1e-3,N,α)
end
function Base.zero(prior::ConstIsoPrior,C::Coherence)
  d = prod(size(C,2,3))
  IsoMultiNormalStats(fill(zero(prior.S),d),fill(prior.S,d),
                      prior.N,prior.N,prior.α)
end

function Base.:(+)(a::ConstIsoPrior{T},b::IsoMultiNormalStats{T}) where T
  d = length(b.μ)
  n1 = b.n1 + a.N
  n2 = b.n2 + a.N
  IsoMultiNormalStats(b.μ.*(b.n1/n1),a.S.+b.S,n1,n2,a.α+b.α)
end

function logpdf(stats::ConstIsoPrior,x::AbstractVector)
  α = stats.N/2 + stats.α
  β = 0.5stats.S
  σ = β./α .* ((stats.N + 1) ./ stats.N)

  sum(logpdf.(TDist(2α),x./σ))
end
