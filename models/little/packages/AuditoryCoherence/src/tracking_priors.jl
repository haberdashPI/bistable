import Distributions: Normal, logpdf

mutable struct MultiNormalStats{T}
  μ::Vector{T}
  S::Matrix{T}
  n::Float64
  x2_offset::Float64
end
function MultiNormalStats(data::AbstractMatrix,n=size(data,1),
                          x2_offset=size(data,2))
  μ = mean(data,1)
  S = (data.-μ)'(data.-μ) .* n./size(data,1)
  MultiNormalStats{eltype(data)}(μ,S,n,x2_offset)
end
MultiNormalStats(::Type{T},n) where T =
  MultiNormalStats{T}(zeros(T,n),zeros(T,n,n),0,0)

function update!(stats::MultiNormalStats{T},x::AbstractVector{T}) where T
  δ = x .- stats.μ
  stats.n += 1
  stats.μ .+= δ/stats.n

  δ2 = x .- stats.μ
  stats.x2 .+= δ.*δ2'

  stats
end
update(stats::MultiNormalStats{T},x::AbstractVector{T}) where T =
  update!(copy(stats),x)

function mult!(stats::MultiNormalStats,c)
  stats.n *= c

  stats
end

function Base.:(+)(a::MultiNormalStats{T},b::MultiNormalStats{T}) where T
  N = a.n + b.n
  aw = a.n/N
  bw = b.n/N
  MultiNormalStats(a.μ.*aw .+ b.μ.*bw,a.S.*aw + b.S.*bw,N,
                   a.x2_offset + b.x2_offset)
end

function logpdf(stats::MultiNormalStats,x::AbstractVector)
  d = length(x)
  Λ = stats.S ./ stats.n
  μ = stats.μ
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

function logpdf_thresh(stats::MultiNormalStats,x::AbstractVector,thresh)
  ii = find(diag(stats.x2) .- (stats.x .* stats.x) .> thresh)

  d = length(ii)
  Λ = stats.S[ii,ii] ./ stats.n
  μ = stats.μ[ii]
  n = stats.n
  # a bit of cheating...
  nu = max(1,floor(Int,stats.n + stats.x2_offset - d + 1))

  jj = setdiff(1:length(x),ii)
  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x[ii]) +
    sum(logpdf.(Normal(0,thresh),x[jj] .- stats.x[jj]))
end

mutable struct IsoMultiNormalStats{T}
  μ::Vector{T}
  S::Vector{T}
  n::Float64
  α::Float64
end

function IsoMultiNormalStats(data::AbstractMatrix,n=size(data,1),α=1,β=1)
  μ = mean(data,1)
  n = size(data,1)
  S = data.^2
  IsoMultiNormalStats{eltype(data)}(μ,S,n,α,β)
end

function update!(stats::IsoMultiNormalStats{T},x::AbstractVector{T}) where T
  δ = x .- stats.μ
  stats.n += 1
  stats.μ .+= δ/stats.n

  δ2 = x .- stats.μ
  stats.S .+= δ.*δ2

  stats
end
update(stats::IsoMultiNormalStats{T},x::AbstractVector{T}) where T =
  update!(copy(stats),x)

function mult!(stats::IsoMultiNormalStats,c)
  stats.n *= c

  stats
end

function Base.:(+)(a::IsoMultiNormalStats{T},b::IsoMultiNormalStats{T}) where T
  N = a.n + b.n
  IsoMultiNormalStats(a.μ.*(a.n./N) .+ b.μ.*(b.n./N),a.S.+b.S,
                      N,a.α+b.α)
end

function logpdf(stats::IsoMultiNormalStats,x::AbstractMatrix)
  μ = stats.μ
  α = stats.n + stats.α
  β = stats.S
  σ = β./α .* ((stats.n + 1) ./ stats.n)

  logpdf(TDist(2α),(x .- μ)./σ)
end

function logpdf_thresh(stats::IsoMultiNormalStats,x::AbstractMatrix,thresh)
  ii = find(diag(stats.x2) .- (stats.x .* stats.x) .> thresh)

  μ = stats.μ[ii]
  α = stats.n + stats.α
  β = stats.S[ii,ii]
  σ = β./α .* ((stats.n + 1) ./ stats.n)

  logpdf(TDist(2α),(x .- μ)./σ)
end

