export ridgenorm, rdist

function rdist(;scale=nothing,freq=nothing)
  (a,b) -> (a[1] - b[1])^2 / (scale^2) + (a[2] - b[2])^2 / (freq^2)
end

mutable struct RidgeMultiNormalStats{T} <: Stats{T}
  μ::Vector{T}
  S::Vector{T}
  n::Float64
  x2_offset::Float64
  corr::SparseMatrixCSC{T,Int}
end
Base.std(x::RidgeMultiNormalStats) = x.n > 0 ? x.S ./ x.n : Inf
Base.zero(x::RidgeMultiNormalStats{T},C::Coherence) where T =
  RidgeMultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0,x.corr)

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

function ridgenorm(prior::Coherence,n::Number,
                   x2_offset::Number=prod(size(prior,2,3));
                   thresh=1e-3,scale=nothing,
                   freq=nothing)
  ridgenorm(prior,n,rdist(scale=scale,freq=freq),x2_offset,thresh=thresh)
end

function ridgenorm(prior::Coherence,n::Number,dist::Function,
                   x2_offset::Number=prod(size(prior,2,3));thresh=1e-3)
  data = reshape(mean(prior,4),size(prior,1),:)
  μ = squeeze(mean(data,1),1)
  S = squeeze(sum(data.^2,1),1) .* n./size(data,1) .+ n.*thresh.^2
  corr = findcorr(size(prior,2,3),dist)
  RidgeMultiNormalStats{eltype(data)}(μ,S,n,x2_offset,corr)
end

# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
function update!(stats::RidgeMultiNormalStats{T},x::AbstractVector{T},
                 w=1.0) where T
  stats.n += w

  μ₀ = copy(stats.μ)
  stats.μ .+= (x .- stats.μ).*(w ./ stats.n)
  stats.S .+= (x .- stats.μ).*(x .- μ₀).*w

  stats
end

function downdate!(stats::RidgeMultiNormalStats{T},x::AbstractVector{T},
                   w=1.0) where T
  if w ≈ stats.n
    stats.μ .= 0
    stats.S .= 0
    stats.n = 0
  elseif w < stats.n
    μₙ = copy(stats.μ)
    stats.μ .= (stats.μ .- x.*(w ./ stats.n)) ./ (1 .- (w ./ stats.n))
    stats.S .-= (x .- stats.μ).*(x .- μₙ).*w

    stats.n -= w
  else
    error("Sum of weights would become negative. This happens when trying ",
          "to remove more samples then were added using `update!`.")
  end

  stats
end

function mult!(stats::RidgeMultiNormalStats,c)
  stats.n *= c

  stats
end

function Base.:(+)(a::RidgeMultiNormalStats{T},b::RidgeMultiNormalStats{T}) where T
  n = a.n + b.n
  RidgeMultiNormalStats(a.μ.*(a.n./n) .+ b.μ.*(b.n./n), a.S.+b.S,
                        n,a.x2_offset+b.x2_offset,a.corr)
end

function logpdf(stats::RidgeMultiNormalStats,x::AbstractVector)
  d = length(stats.μ)
  σ = sqrt.(stats.S ./ stats.n)
  Λ = (Diagonal(σ) * stats.corr) * Diagonal(σ)
  nu = stats.n + stats.x2_offset - d + 1

  if nu < 1
    error("Insufficient data for predictive distribution, ",
          "try collecting more data or increasing x2_offset")
  end

  logpdf_mvt(nu,stats.μ,Λ .* ((stats.n + 1)/(stats.n * nu)),x)
end

struct ConstRidgePrior{T} <: Stats{T}
  S::T
  N::Float64
  x2_offset::Float64
  corr::SparseMatrixCSC{T,Int}
end

function ridgenorm(prior::Number,N::Number,dims::Tuple,dist::Function,
                   x2_offset::Number=prod(dims))
  ConstRidgePrior{typeof(prior)}(prior,N,x2_offset,findcorr(dims,dist))
end

function ridgenorm(prior::Number,N::Number,dims::Tuple,
                   x2_offset::Number=prod(dims);
                   scale=nothing,freq=nothing)
  dist = rdist(scale=scale,freq=freq)
  ridgenorm(prior,N,dims,dist,x2_offset)
end

function Base.zero(prior::ConstRidgePrior,C::Coherence)
  d = prod(size(C,2,3))
  RidgeMultiNormalStats(fill(zero(prior.S),d),fill(prior.S,d),
                        prior.N,prior.x2_offset,prior.corr)
end

function Base.:(+)(a::ConstRidgePrior{T},b::RidgeMultiNormalStats{T}) where T
  d = length(b.μ)
  n = b.n + a.N
  RidgeMultiNormalStats(b.μ.*(b.n/n),a.S.+b.S,n,a.x2_offset+b.x2_offset,
                        a.corr)
end

function logpdf(stats::ConstRidgePrior,x::AbstractVector)
  d = size(stats.S,1)
  σ = stats.S / stats.N
  Λ = stats.corr .* σ
  nu = max(1,floor(Int,stats.N + stats.x2_offset - d + 1))

  logpdf_mvt(nu,0,Λ .* ((stats.N + 1)/(stats.N * nu)),x)
end
