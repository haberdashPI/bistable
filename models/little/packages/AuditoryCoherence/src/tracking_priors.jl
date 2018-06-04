import Distributions: TDist, Normal, logpdf

struct Beta
  α::Float64
  β::Float64
end
update(d::Beta,x::Bool,w=1.0) = x ? Beta(d.α+w,d.β) : Beta(d.α,d.β+w)
mult(d::Beta,x) = Beta(d.α*x,d.β*x)
Base.:(+)(a::Beta,b::Beta) = Beta(a.α+b.α,a.β+b.β)

struct BinomialCond{T}
  at::Dict{T,Beta}
end
BinomialCond(pairs::Pair{T}...) where T = BinomialCond{T}(Dict(pairs...))
update!(d::BinomialCond,cond,x::Bool,w=1.0) = d.at[cond] = update(d.at[cond],x,w)
function mult!(d::BinomialCond,x)
  for cond in keys(d.at)
    d.at[cond] = mult(d.at[cond],x)
  end
end
Base.:(+)(a::BinomialCond,b::BinomialCond) =
  BinomialCond(merge(+,a.at,b.at))
function logpdf(d::BinomialCond,cond,x)
  p = d.at[cond].α / (d.at[cond].α + d.at[cond].β)
  x ? log(p) : log(1-p)
end

mutable struct MultiNormalStats{T}
  μ::Vector{T}
  S::Matrix{T}
  n1::Float64
  n2::Float64
  x2_offset::Float64
end
Base.zero(x::MultiNormalStats{T}) where T =
  MultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0)
function MultiNormalStats(data::AbstractMatrix,n=size(data,1),
                          x2_offset=size(data,2))
  μ = mean(data,1)
  S = (data.-μ)'(data.-μ) .* n./size(data,1)
  MultiNormalStats{eltype(data)}(μ,S,n,n,x2_offset)
end
MultiNormalStats(::Type{T},n) where T =
  MultiNormalStats{T}(zeros(T,n),zeros(T,n,n),0,0)

function update!(stats::MultiNormalStats{T},x::AbstractVector{T},w=1.0) where T
  δ = x .- stats.μ
  stats.n1 += w
  stats.μ .+= δ * (w/stats.n2)

  δ2 = x .- stats.μ
  stats.n2 += w^2
  stats.S .+= δ.*δ2' * (w2/stats.n2)

  stats
end
update(stats::MultiNormalStats{T},x::AbstractVector{T}) where T =
  update!(copy(stats),x)

function mult!(stats::MultiNormalStats,c)
  stats.n1 *= c
  stats.n2 *= c^2

  stats
end

function Base.:(+)(a::MultiNormalStats{T},b::MultiNormalStats{T}) where T
  N1 = a.n1 + b.n1
  aw = a.n1/N1
  bw = b.n2/N1
  μ = a.μ.*aw .+ b.μ.*bw

  N2 = a.n2 + b.n2
  aw = a.n2/N2
  bw = b.n2/N2
  S = a.S.*aw + b.S.*bw
  MultiNormalStats(μ,S,N1,N2,a.x2_offset + b.x2_offset)
end

function logpdf(stats::MultiNormalStats,x::AbstractVector)
  d = length(x)
  Λ = stats.S ./ stats.n2
  μ = stats.μ
  n = stats.n1
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n2 + stats.x2_offset - d + 1))

  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x)
end
pdf(stats::MultiNormalStats,x::AbstractVector) = exp(logpdf(stats,x))

function logpdf_mvt(v,μ,Σ,x)
  d = length(μ)
  C = lgamma((v+1)/2) - (lgamma(v/2)+log(v*π)^(d/2)) - 0.5logabsdet(Σ)[1]
  diff = abs.(μ.-x)

  C*log(1+1/v*(diff'/Σ)*diff)*-(v+d)/2
end

function logpdf_thresh(stats::MultiNormalStats,x::AbstractVector,thresh)
  ii = find(stats.S ./ stats.n .> thresh)

  d = length(ii)
  Λ = stats.S[ii,ii] ./ stats.n2
  μ = stats.μ[ii]
  n = stats.n1
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n2 + stats.x2_offset - d + 1))

  jj = setdiff(1:length(x),ii)
  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x[ii]) +
    sum(logpdf.(Normal(0,thresh),x[jj] .- stats.μ[jj]))
end

mutable struct IsoMultiNormalStats{T}
  μ::Vector{T}
  S::Vector{T}
  n1::Float64
  n2::Float64
  α::Float64
end
Base.std(x::IsoMultiNormalStats) = x.S ./ x.n2
Base.zero(x::IsoMultiNormalStats{T}) where T =
  IsoMultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0,0.0)

function IsoMultiNormalStats(data::AbstractMatrix,n=size(data,1),α=1)
  μ = squeeze(mean(data,1),1)
  S = squeeze(sum(data.^2,1),1) .* n./size(data,1)
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
update(stats::IsoMultiNormalStats{T},x::AbstractVector{T}) where T =
  update!(copy(stats),x)

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

function logpdf_thresh(stats::IsoMultiNormalStats,x::AbstractVector,thresh)
  ii = find(stats.S ./ stats.n2 .> thresh)

  μ = stats.μ[ii]
  α = stats.n1/2 + stats.α
  β = 0.5stats.S[ii]
  σ = β./α .* ((stats.n2 + 1) ./ stats.n2)

  jj = setdiff(1:length(x),ii)
  sum(logpdf.(TDist(2α),(x[ii] .- μ)./σ)) +
    sum(logpdf.(Normal(0,thresh),x[jj] .- stats.μ[jj]))
end

