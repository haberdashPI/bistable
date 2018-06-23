export ridgenorm, rdist

function rdist(;scale=nothing,freq=nothing)
  (a,b) -> (a[1] - b[1])^2 / (scale^2) + (a[2] - b[2])^2 / (freq^2)
end

mutable struct RidgeMultiNormalStats{T} <: Stats{T}
  μ::Vector{T}
  S::Vector{T}
  n1::Float64
  n2::Float64
  x2_offset::Float64
  corr::SparseMatrixCSC{T,Int}
end
Base.std(x::RidgeMultiNormalStats) = x.n2 > 0 ? x.S ./ x.n2 : Inf
Base.zero(x::RidgeMultiNormalStats{T},C::Coherence) where T =
  RidgeMultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0,0.0,x.corr)

function findcorr(dims,dist)
  n = prod(dims)
  corr = spzeros(n,n)
  for (i,ii) in enumerate(CartesianRange(dims))
    for (j,jj) in enumerate(CartesianRange(dims))
      d = dist(ii.I,jj.I)
      if d < 3^2
        corr[i,j] = exp(-d^2)
      end
    end
  end

  corr
end

function ridgenorm(prior::Coherence,n,x2_offset=1;scale=nothing,freq=nothing)
  ridgenorm(prior,n,rdist(;scale=scale,freq=freq),x2_offset)
end

function ridgenorm(prior::Coherence,n,dist,x2_offset =1)
  data = reshape(mean(prior,4),size(prior,1),:)
  μ = squeeze(mean(data,1),1)
  S = squeeze(sum(data.^2,1),1) .* n./size(data,1)
  corr = findcorr(size(prior,2,3),dist)
  RidgeMultiNormalStats{eltype(data)}(μ,S,n,n,x2_offset,corr)
end

function update!(stats::RidgeMultiNormalStats{T},x::AbstractVector{T},w=1.0) where T
  δ = x .- stats.μ
  stats.n1 += w
  stats.μ .+= δ * (w/stats.n1)

  δ2 = x .- stats.μ
  stats.n2 += w^2
  stats.S .+= δ.*δ2 * (w^2/stats.n2)

  stats
end

function mult!(stats::RidgeMultiNormalStats,c)
  stats.n1 *= c
  stats.n2 *= c^2

  stats
end

function Base.:(+)(a::RidgeMultiNormalStats{T},b::RidgeMultiNormalStats{T}) where T
  n1 = a.n1 + b.n1
  n2 = a.n2 + b.n2
  RidgeMultiNormalStats(a.μ.*(a.n1./n1) .+ b.μ.*(b.n1./n1), a.S.+b.S,
                        n1,n2,a.x2_offset+b.x2_offset,a.corr)
end


function logpdf_thresh(stats::RidgeMultiNormalStats,x::AbstractVector,thresh)
  ii = find(stats.S ./ stats.n2 .> thresh)

  d = length(ii)
  μ = stats.μ[ii]
  σ = sqrt.(stats.S[ii] ./ stats.n2)
  Λ = (Diagonal(σ) * stats.corr[ii,ii]) * Diagonal(σ)
  n = stats.n1
  nu = max(1,floor(Int,stats.n2 + stats.x2_offset - d + 1))

  jj = setdiff(1:length(x),ii)
  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x[ii]) +
    sum(logpdf.(Normal(0,2thresh),x[jj] .- stats.μ[jj]))
end

struct ConstRidgePrior{T} <: Stats{T}
  S::T
  N::Float64
  x2_offset::Float64
  corr::SparseMatrixCSC{T,Int}
end
function ridgenorm(prior::Number,N,dims,dist,x2_offset=1)
  ConstRidgePrior{typeof(prior)}(prior,N,x2_offset,findcorr(dims,dist))
end
function Base.zero(prior::ConstRidgePrior,C::Coherence)
  RidgeMultiNormalStats(fill(zero(prior.S),d),fill(prior.S,d),
                        prior.N,prior.N,prior.x2_offset,
                        prior.corr)
end

function Base.:(+)(a::ConstRidgePrior{T},b::RidgeMultiNormalStats{T}) where T
  d = length(b.μ)
  n1 = b.n1 + a.N
  n2 = b.n2 + a.N
  RidgeMultiNormalStats(b.μ.*(b.n1/n1),a.S.+b.S,n1,n2,a.x2_offset+b.x2_offset,
                       prior.corr)
end

function logpdf_thresh(stats::ConstRidgePrior,x::AbstractVector,thresh)
  μ = 0
  σ = sqrt(stats.S / stats.N)
  Λ = (Diagonal(σ) * stats.corr) * Diagonal(σ)
  n = stats.N
  nu = max(1,floor(Int,stats.N + stats.x2_offset - d + 1))

  logpdf_mvt(nu,μ,Λ .* ((n + 1)/(n * nu)),x)
end
