# TODO: make eigenspace time unit aware and rename it
struct EigenSpace{U,L}
  u::Matrix{U}
  λ::Vector{L}
  var::Vector{L}
end
ncomponents(x::EigenSpace) = size(x.u,2)

function EigenSpace(u::Array,λ::Array{L},var::Array{L}) where L
  @assert size(u,2) == size(λ,1) "Mismatched number of eigenvectors and values."
  EigenSpace(u[:,:],λ[:],var[:])
end

EigenSpace(d::Number,n::Number) = EigenSpace(Float64,d,n)
EigenSpace(::Type{T},d::Number,n::Number) where T =
  EigenSpace([I; fill(zero(T),d-n,n)],fill(one(real(T)),n),ones(real(T),d))

Base.show(io::IO,x::EigenSpace) = write(io,"EigenSpace(λ = ",string(x.λ),")")

# TODO: use last index as time dimension
struct EigenSeries{U,L} <: AbstractArray{EigenSpace{U,L},1}
  u::Array{U,3}
  λ::Array{L,2}
  var::Array{L,2}
  delta::Seconds{Float64}
end
ncomponents(x::EigenSeries) = size(x.u,3)

EigenSeries(::Type{T},t,d,n,delta) where T =
  EigenSeries(Array{T}(t,d,n),Array{real(T)}(t,n),Array{real(T)}(t,d),delta)
EigenSeries(t,d,n,delta) = EigenSeries(Float64,t,d,n,delta)
EigenSeries(t,x::EigenSpace) =
  EigenSeries(similar(x.u,(t,size(x.u)...)),
              similar(x.λ,(t,size(x.λ)...)),
              similar(x.λ,(t,size(x.u,1))),
              x.delta)

function Base.similar(x::EigenSeries,::Type{EigenSpace{U,L}},
                      dims::NTuple{N,Int}) where {U,L,N}

  @assert length(dims) == 1 "EigenSeries can only have a dimensionality of 1."
  EigenSeries(similar(x.u,U,(dims[1],size(x.u,2),size(x.u,3))),
              similar(x.λ,L,(dims[1],size(x.λ,2))),
              similar(x.var,L,(dims[1],size(x.u,1))),
              x.delta)
end

function Base.setindex!(x::EigenSeries,v::EigenSpace,i::Int)
  x.u[i,:,:] = v.u
  x.λ[i,:] = v.λ
  x.var[i,:] = v.var
  v
end
EigenSpace(x::EigenSeries) = EigenSpace(size(x.u,2),size(x.u,3))
Base.getindex(x::EigenSeries,i::Int) = EigenSpace(x.u[i,:,:],x.λ[i,:],x.var[i,:])
Base.getindex(x::EigenSeries,i::Quantity{N,TimeDim}) where N =
  x[max(1,floor(Int,i / x.delta))]
Base.size(x::EigenSeries) = (size(x.u,1),)
Base.IndexStyle(::EigenSeries) = IndexLinear()

update_constants(n::Int) = n/(n+1),1/(n+1)
update_constants(f::Float64) = (1-f),f
update(pca::EigenSpace,x::AbstractArray,c) = update(pca,vec(x),c)

function update(pca::EigenSpace,x::AbstractVector,c)
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

  var = c_o.*var + c_n.*abs2.(x)
  EigenSpace(u[:,indices],λ[indices],var)
end
Base.eigvals(pca::EigenSpace) = pca.λ
Base.eigvecs(pca::EigenSpace) = pca.u
Base.var(pca::EigenSpace) = pca.var

# aproximate y in terms of eigenvectors of x
function project(x::EigenSpace,y::EigenSpace)
  if x.u === y.u
    y
  else
    λp = sum((x.u'y.u).^2 .* y.λ',2)
    EigenSpace(x.u,λp,x.n)
  end
end
