using Unitful
using NMF

struct NMFSpace <: Factors
  W::Matrix{Float64}
  H::Matrix{Float64}
  delta::Seconds{Float64}

  function NMFSpace(W::Matrix{Float64},H::Matrix{Float64},
                    delta::Seconds{Float64})
    @assert size(W,2) == size(H,1)
    new(W,H,delta)
  end
end
nslices(x::NMFSpace) = size(x.W,1)
ncomponents(x::NMFSpace) = size(x.W,2)
nunits(x::NMFSpace) = size(x.H,2)

NMFSpace(t::Number,d::Number,n::Number) = NMFSpace(zeros(t,n),zeros(n,d))
Base.show(io::IO,x::NMFSpace) = write(io,"NMFSpace(k=$(ncomponents(x)))")

factors(x::NMFSpace) = x.H'
strengths(x::NMFSpace) = mean(x.W,1)

struct NMFSeries <: FactorSeries{NMFSpace}
  W::Array{Float64,3}
  H::Array{Float64,3}
  delta::Seconds{Float64}
  step::Seconds{Float64}

  function NMFSeries(W::Array{Float64,3},H::Array{Float64,3},
                     delta::Seconds{Float64},step::Seconds{Float64})
    @assert size(W,3) == size(H,2)
    @assert size(W,1) == size(H,1)
    @assert delta <= step
    new(W,H,delta,step)
  end
end
ncomponents(x::NMFSeries) = size(x.W,3)
nunits(x::NMFSeries) = size(x.H,3)
nslices(x::NMFSeries) = size(x.W,2)

NMFSeries(t,w,d,k,delta,step) =
  NMFSeries(zeros(t,w,k),zeros(t,k,d),delta,step)
NMFSeries(t,step,x::NMFSpace) =
  NMFSeries(t,nslices(x),nunits(x),ncomponents(x),x.delta,step)

function Base.similar(x::NMFSeries,::Type{NMFSeries},
                      dims::NTuple{N,Int}) where {N}
  @assert length(dims) == 1 "NMFSeries can only have dimensionality of 1."
  NMFSeries(similar(x.W,(dims[1],size(x.W,2),size(x.W,3))),
            similar(x.H,(dims[1],size(x.H,2),size(x.H,3))),
            x.delta,x.step)
end

function Base.setindex!(x::NMFSeries,v::NMFSpace,i::Int)
  @assert nslices(x) >= nslices(v)
  @assert nunits(x) == nunits(v)
  @assert ncomponents(x) == ncomponents(v)

  x.W[i,1:size(v.W,1),:] = v.W
  x.W[i,size(v.W,1)+1:end,:] = 0
  x.H[i,:,:] = v.H
  v
end
NMFSpace(x::NMFSeries) = NMFSpace(nslices(x),ncomponents(x),nunits(x))
Base.getindex(x::NMFSeries,i::Int) = NMFSpace(x.W[i,:,:],x.H[i,:,:],x.delta)
Base.getindex(x::NMFSeries,i::Quantity{N,TimeDim}) where N =
  x[max(1,floor(Int,i / x.step))]
Base.size(x::NMFSeries) = (size(x.W,1),)
Base.IndexStyle(::NMFSeries) = IndexLinear()

factors(x::NMFSeries) = permutedims(x.H,[1,3,2])
strengths(x::NMFSeries) = vec(mean(x.W,(1,2))) .* vec(mean(x.H,(1,3)))

function bestordering(x,y)
  @assert length(x) == length(y)
  selected = collect(1:length(y))
  for i in 1:length(y)-1
    k = indmin(norm(x[i] .- y[j]) for j in selected[i:end])
    selected[i],selected[i+k-1] = selected[i+k-1],selected[i]
  end

  selected
end


function normalize_components(C::NMFSeries)
  K = ncomponents(C)

  # rearrange to maintain index for similar components
  # TODO: take some short leaky average of the components
  Ĥ = C.H[1,:,:]
  α = 0.1
  for t in 2:size(C.H,1)
    ordering = bestordering([Ĥ[k,:] for k in 1:K],
                            [C.H[t,k,:] for k in 1:K])

    C.H[t,:,:] = C.H[t,ordering,:]
    C.W[t,:,:] = C.W[t,:,ordering]
    Ĥ .+= (1-α)Ĥ + α*C.H[t,:,:]
  end

  # rearrange so largest is first
  ordering = sortperm(strengths(C),rev=true)
  C.H .= C.H[:,ordering,:]
  C.W .= C.W[:,:,ordering]

  C
end
