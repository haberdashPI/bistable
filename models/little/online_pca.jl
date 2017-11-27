struct CCI_PCA{T}
  v::Matrix{T}
  n::Int
  ε::Float64
end

struct IPCA{T}
  u::Matrix{T}
  λ::Vector{T}
  n::Int
end

const default_pca_method = :ipca

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

function Base.push!(pca::CCI_PCA,x::Vector)
  v,n,ε = pca.v,pca.n,pca.ε

  total = norm(x)
  x = copy(x)
  for i in 1:size(v,2)
    if norm(x) < pca.ε * total
      break
    end

    λ_i = norm(v[:,i])
    if iszero(λ_i)
      v[:,i] = x
      break
    else
      u_i = v[:,i] ./ λ_i
      xu_i = dot(x,u_i)
      v[:,i] .= n/(n+1).*v[:,i] .+ (1/(n+1)).*(x .* xu_i)
      x .-= xu_i .* u_i
    end
  end

  CCI_PCA(v,n+1,pca.ε)
end

normalize_(x) = (n = norm(x); iszero(n) ? x : x./n)
Base.eigvals(pca::CCI_PCA) = vec(mapslices(norm,pca.v,1))
Base.eigvecs(pca::CCI_PCA) = mapslices(normalize_,pca.v,1)

function Base.push!(pca::IPCA,x::Vector)
  u,λ,n = pca.u,pca.λ,pca.n
  d = length(λ)

  x_u = pca.u'x
  x_p = x - pca.u*x_u
  nx_p = norm(x_p)

  Q = similar(u,d+1,d+1)
  Q[1:d,1:d] .= (n+1).*Diagonal(λ) .+ x_u.*x_u'
  Q[1:d,d+1] .= Q[d+1,1:d] .= nx_p.*x_u
  Q[d+1,d+1] = nx_p^2
  Q .*= n/(n+1)^2

  λ,v = eig(Q)
  u = [u x_p./nx_p] * v

  _,mini = findmin(abs.(λ))
  indices = vcat(1:mini-1,mini+1:d+1)

  IPCA(u[:,indices],λ[indices],n+1)
end
Base.eigvals(pca::IPCA) = pca.λ
Base.eigvecs(pca::IPCA) = pca.u
