using Einsum

abstract type OnlinePCA{T} end

struct CCI_PCA{T} <: OnlinePCA{T}
  v::Matrix{T}
  n::Int
  ε::Float64
end
ncomponents(x::CCI_PCA) = size(v,2)

struct IPCA{T} <: OnlinePCA{T}
  u::Matrix{T}
  λ::Vector{T}
  n::Int
end
ncomponents(x::CCI_PCA) = size(u,2)

const default_pca_method = :ccipca

# initialize using dimensions and identity matrix
OnlinePCA(d::Number,n::Number;ε=1e-2,method=default_pca_method) =
  method == :ccipca ? CCI_PCA([I; fill(0.0,d-n,n)],1,ε) :
  method == :ipca ? IPCA([I; fill(0.0,d-n,n)],fill(1.0,n),1) :
  error("No online-pca method named \"$method\".")

# initialize using singular values
OnlinePCA(sv::LinAlg.SVD,n::Number;ε=1e-2,method=default_pca_method) =
  method == :ccipca ? CCI_PCA(sv[:V] * Diagonal(sv[:S].^2),n,ε) :
  method == :ipca ? IPCA(sv[:V],sv[:S].^2,n) :
  error("No online-pca method named \"$method\".")

update_constants(::Void,n) = n/(n+1),1/(n+1)
update_constants(f::Number,n) = (1-f),f

function Base.LinAlg.lowrankupdate(pca::CCI_PCA,x::AbstractVector,f=nothing)
  v,n,ε = copy(pca.v),pca.n,pca.ε

  total = norm(x)
  x = copy(x)
  for i in 1:size(v,2)
    if norm(x) < pca.ε * total
      break
    end

    λ_i = norm(v[:,i])
    c_o,c_n = update_constants(n,f)
    if iszero(λ_i)
      v[:,i] = x
      break
    else
      u_i = v[:,i] ./ λ_i
      xu_i = dot(x,u_i)
      v[:,i] .= c_o.*v[:,i] .+ c_n.*(x .* xu_i)
      x .-= xu_i .* u_i
    end
  end

  CCI_PCA(v,n+1,pca.ε)
end

normalize_(x) = (n = norm(x); iszero(n) ? x : x./n)
Base.eigvals(pca::CCI_PCA) = vec(mapslices(norm,pca.v,1))
Base.eigvecs(pca::CCI_PCA) = mapslices(normalize_,pca.v,1)

function Base.LinAlg.lowrankupdate(pca::IPCA,x::AbstractVector,f=nothing)
  u,λ,n = pca.u,pca.λ,pca.n
  d = length(λ)

  x_u = u'x
  x_p = x - u*x_u
  nx_p = norm(x_p)

  c_o,c_n = update_constants(f,n)
  Q = similar(u,d+1,d+1)
  Q[1:d,1:d] .= c_o.*(Diagonal(λ) .+ c_n.*x_u.*x_u')
  Q[1:d,d+1] .= Q[d+1,1:d] .= c_o*c_n .* nx_p.*x_u
  Q[d+1,d+1] = c_o*c_n * nx_p^2

  λ,v = eig(Q)
  u = [u x_p./nx_p] * v

  min_i = indmin(abs(λ[i]) for i in eachindex(λ))
  indices = vcat(1:min_i-1,min_i+1:d+1)

  IPCA(u[:,indices],λ[indices],n+1)
end
Base.eigvals(pca::IPCA) = pca.λ
Base.eigvecs(pca::IPCA) = pca.u

Base.:(-)(a::IPCA,b::IPCA) = op_helper(-,a,b)
Base.:(+)(a::IPCA,b::IPCA) = op_helper(+,a,b)
function op_helper(op,x::IPCA,y::IPCA)
  IPCA(x.u,op(x.λ,project(x,y)),x.n)
end

# aproximate eigenvalues of y
# in terms of eigenvalues of x
function project(x::IPCA,y::IPCA)
  a = x.u'y.u

  # λp = diag(a * Diagonal(y.λ) * a')
  @einsum λp[i] := a[i,j]^2*x.λ[j]
  λp
end

# NOTE: to get this
# working I need the second half of the
# algorithm to allow for reduced
# rank approximations...
# probably not worth implementing
#========================================
function gram_schmit(U,X)
  Up = similar(U,size(U,1),size(U,2) + size(X,2))
  n = size(U,2)

  projn(u,v) = dot(u,v)*u
  function projsum(i)
    sum(1:n+i-1) do j
      if j <= n
        projn(U[:,j],X[:,i])
      else
        projn(Xp[:,j-n],X[:,i])
      end
    end
  end

  for i in size(X,2)
    Xp[:,i] = normalize(X[:,i] .- projsum(i))
  end

  Xp
end

function blockupdate(pca::RSVD,X::AbstractMatrix{T}) where T
  U,S,V = pca.U,pca.S,pca.V

  # U update
  m = size(U,2)
  Up = gram_schmidt(U,X)
  Xp = Up[:,m+1:end]

  # S update
  n,m = size(S)
  Sp = fill(zero(T),size(S) .+ size(X,2)...)
  Sp[1:n,1:m] = S
  Sp[1:n,m+1:end] = U'X
  Sp[n+1:end,m+1:end] = Xp'X

  # V update
  n,m = size(V)
  Vp = fill(zero(T),size(V) .+ size(X)...)
  Vp[1:n,1:m] = V
  Vp[n+1:end,m+1:end] = I

  # recaluate SVD
  sv, = svds(Sp,nsv=pca.ncomponents)

  RSVD(Up*sv[:U],sv[:S],sv[:Vt]*Vp')
end
========================================#
