using Einsum

# TODO??: implement EigenSeries to work like EigenSpace, but with the memory managed a
# little more efficiently
struct EigenSeries{T}
  u::Array{T,3}
  λ::Array{T,2}
end
EigenSeries(x,n) =
  EigenSeries(similar(x,size(x,1),prod(size(x)[2:end]),n),
              similar(x,size(x,1),n))

struct EigenSpace{T}
  u::Matrix{T}
  λ::Vector{T}
end
ncomponents(x::EigenSpace) = size(x.u,2)
EigenSpace(x::EigenSeries{T}) EigenSpace(size(x.u,2),size(x.u,3))

Base.setindex!(x::EigenSeries,v::EigenSpace,i::Int) =
  EigenSpace(x.u[i,:,:],x.v[i,:]); v

# initialize using dimensions and identity matrix
EigenSpace(d,n) = EigenSpace([I; fill(0.0,d-n,n)],fill(1.0,n),1)

update_constants(n::Int) = n/(n+1),1/(n+1)
update_constants(f::Float64) (1-f),f

function update(pca::EigenSpace,x::AbstractVector,c)
  u,λ,n = pca.u,pca.λ
  d = length(λ)

  x_u = u'x
  x_p = x - u*x_u
  nx_p = norm(x_p)

  c_o,c_n = update_constants(c)
  Q = similar(u,d+1,d+1)
  Q[1:d,1:d] .= c_o.*(Diagonal(λ) .+ c_n.*x_u.*x_u')
  Q[1:d,d+1] .= Q[d+1,1:d] .= c_o*c_n .* nx_p.*x_u
  Q[d+1,d+1] = c_o*c_n * nx_p^2

  λ,v = eig(Q)
  u = [u x_p./nx_p] * v

  min_i = indmin(abs(λ[i]) for i in eachindex(λ))
  indices = vcat(1:min_i-1,min_i+1:d+1)

  EigenSpace(u[:,indices],λ[indices],n+1)
end
Base.eigvals(pca::EigenSpace) = pca.λ
Base.eigvecs(pca::EigenSpace) = pca.u

Base:(-)(a::EigenSpace,b::EigenSpace) = op_helper(-,a,b)
addp(a::EigenSpace,b::EigenSpace) = op_helper(+,a,b)
function op_helper(op,x::EigenSpace,y::EigenSpace)
  EigenSpace(x.u,op(x.λ,project(x,y)),x.n)
end

# aproximate y in terms of eigenvectors of x
function project!(x::EigenSpace,y::EigenSpace)
  if x.u === y.u
    y
  else
    p = x.u'y.u

    # λp = diag(p * Diagonal(y.λ) * p')
    @einsum λp[i] := p[i,j]^2*y.λ[j]
    EigenSpace(x.u,λp,x.n)
  end
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
