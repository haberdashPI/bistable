abstract type EigenSpace{U,L} end

struct SimpleEigenSpace{U,L} <: EigenSpace{U,L}
  u::Vector{U}
  λ::L
  var::L
end
ncomponents(x::SimpleEigenSpace) = 1
Base.:(*)(x::EigenSpace,y::Array) = (x.u'y.*x.λ)' * x.u'
function Base.convert(::Type{SimpleEigenSpace{U,L}},
                      y::SimpleEigenSpace{U2,L2}) where {U,L,U2,L2}
  SimpleEigenSpace{U,L}(y.u,y.λ,y.var)
end

struct FullEigenSpace{U,L} <: EigenSpace{U,L}
  u::Matrix{U}
  λ::Vector{L}
  var::L
end
ncomponents(x::FullEigenSpace) = size(x.u,2)

function EigenSpace(u::Array,λ::Array{L},var::L) where L
  @assert size(u,2) == size(λ,1) "Mismatched number of eigenvectors and values."
  if size(u,2) == 1
    SimpleEigenSpace(u[:,1],λ[1],var)
  else
    FullEigenSpace(u[:,:],λ[:],var)
  end
end

EigenSpace(d::Number,n::Number) = EigenSpace(Float64,d,n)
EigenSpace(::Type{T},d::Number,n::Number) where T =
  n == 1 ? SimpleEigenSpace(fill(zero(T),d),one(real(T)),one(real(T))) :
  FullEigenSpace([I; fill(zero(T),d-n,n)],fill(one(real(T)),n),one(real(T)))

Base.show(io::IO,x::EigenSpace) = write(io,"EigenSpace(λ = ",string(x.λ),")")

# TODO: use last index as time dimension
struct EigenSeries{U,L} <: AbstractArray{EigenSpace{U,L},1}
  u::Array{U,3}
  λ::Array{L,2}
  var::Array{L,1}
end
ncomponents(x::EigenSeries) = size(x.u,3)

EigenSeries(::Type{T},t,d,n) where T =
  EigenSeries(Array{T}(t,d,n),Array{real(T)}(t,n),Array{real(T)}(t))
EigenSeries(t,d,n) = EigenSeries(Float64,t,d,n)
EigenSeries(t,x::FullEigenSpace) =
  EigenSeries(similar(x.u,(t,size(x.u)...)),
              similar(x.λ,(t,size(x.λ)...)),
              similar(x.λ,(t,)))
EigenSeries(t,x::SimpleEigenSpace) =
  EigenSeries(similar(x.u,(t,length(x.u),1)),
              similar(x.λ,(t,1)),
              similar(x.λ,(t,)))

function Base.similar(x::EigenSeries,::Type{EigenSpace{U,L}},
                      dims::NTuple{N,Int}) where {U,L,N}

  @assert length(dims) == 1 "EigenSeries can only have a dimensionality of 1."
  EigenSeries(similar(x.u,U,(dims[1],size(x.u,2),size(x.u,3))),
              similar(x.λ,L,(dims[1],size(x.λ,2))),
              similar(x.var,L,(dims[1],)))
end

function Base.setindex!(x::EigenSeries,v::EigenSpace,i::Int)
  x.u[i,:,:] = v.u
  x.λ[i,:] = v.λ
  x.var[i] = v.var
  v
end
EigenSpace(x::EigenSeries) = EigenSpace(size(x.u,2),size(x.u,3))
Base.getindex(x::EigenSeries,i::Int) = EigenSpace(x.u[i,:,:],x.λ[i,:],x.var[i])
Base.size(x::EigenSeries) = (size(x.u,1),)
Base.IndexStyle(::EigenSeries) = IndexLinear()

update_constants(n::Int) = n/(n+1),1/(n+1)
update_constants(f::Float64) = (1-f),f
update(pca::EigenSpace,x::AbstractArray,c) = update(pca,vec(x),c)

# TODO: if used, implement residual
function update(pca::SimpleEigenSpace{T,L},x::AbstractVector,c) where {T,L}
  u,λ = pca.u,pca.λ
  c_o,c_n = update_constants(c)
  v = iszero(λ) ? x : c_o.*(λ.*u) .+ c_n.*(x .* x'u)

  λ = norm(v)
  SimpleEigenSpace(vec(v)./λ,λ,zero(L))
end
Base.eigvals(pca::SimpleEigenSpace) = [pca.λ]
Base.eigvecs(pca::SimpleEigenSpace) = pca.u[:,:]
Base.var(pca::SimpleEigenSpace) = pca.var

function update(pca::FullEigenSpace,x::AbstractVector,c)
  u,λ,var = pca.u,pca.λ,pca.var
  c_o,c_n = update_constants(c)

  x_u = u'x
  x_p = x - u*x_u
  nx_p = norm(x_p)

  d = length(λ)
  Q = similar(u,d+1,d+1)
  Q[1:d,1:d] .= c_o.*(Diagonal(λ) .+ c_n.*x_u.*x_u')
  Q[1:d,d+1] .= Q[d+1,1:d] .= c_o*c_n .* nx_p.*x_u
  Q[d+1,d+1] = c_o*c_n * nx_p^2

  λc,v = eig(Q)
  @assert all(iszero.(imag.(λ))) "Unexpected complex eigenvalues"
  λ = real.(λc)
  u = [u x_p./nx_p] * v

  min_i = indmin(abs.(λ))
  indices = vcat(1:min_i-1,min_i+1:d+1)

  var = c_o*var + c_n*sum(abs2.(x))
  FullEigenSpace(u[:,indices],λ[indices],var)
end
Base.eigvals(pca::EigenSpace) = pca.λ
Base.eigvecs(pca::EigenSpace) = pca.u
Base.var(pca::EigenSpace) = pca.var

# aproximate y in terms of eigenvectors of x
function project(x::FullEigenSpace,y::FullEigenSpace)
  if x.u === y.u
    y
  else
    λp = sum((x.u'y.u).^2 .* y.λ',2)
    EigenSpace(x.u,λp,x.n)
  end
end

# TODO: if we go with the SimpleEigenSpace,
# is there some way to think about adaptmi???
