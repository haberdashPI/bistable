using Parameters
using ProgressMeter
include("cortical.jl")

@with_kw struct TCAnalysis <: CorticalRepresentation
  cort::CorticalModel
  τ::Seconds{Float64} = 1s
  sparsity::Float64 = 0.8
  frame_len = 10
  N::Int = 25
  prior::Float64 = 1.0
end
Δt(tc::TCAnalysis) = Δt(tc.cort) * tc.frame_len

numobjects(tc::TCAnalysis) = tc.N
times(tc::TCAnalysis,y) = times(tc.cort,y) .* tc.frame_len

# effecient algorithm to find the quantile of a sparse matrix
# when we know there are no negative entries
function aquantile(x::Union{SparseMatrixCSC,SparseVector},p)
  @assert 0 <= p <= 1 "p must be in the interval [0,1]"
  nz = copy(x.nzval)
  sort!(nz)
  pastzero = false
  inzero = false
  nzeros = length(x) - length(nz)

  real_i = p*length(x)
  start,stop = floor(Int,real_i),ceil(Int,real_i)
  at(i) = i <= nzeros ? zero(eltype(x)) : nz[i - nzeros]
  if start == stop
    at(start)
  elseif stop < length(x)
    interp = real_i - start
    interp*at(start) + (1-interp)*at(stop)
  else
    at(start)
  end
end

function sparsify(x::AbstractVector,p)
  xa = abs.(x)
  cut = quantile(xa,p)
  indices = find(x -> x > cut,xa)
  SparseVector(length(x),indices,x[indices])
end

function sparsify!(C::SparseMatrixCSC,p)
  aC = abs.(C)
  cut = aquantile(aC,p)
  iis = find(0 .< aC.nzval .< cut)
  if length(iis) > 0
    C.nzval[iis] = 0
    dropzeros!(C)
  end
  C
end

matchorder(x,::UniformScaling) = x
function matchorder(x,target)
  matches = x' * target
  order = Array{Int}(size(x,2))
  remaining = IntSet(1:size(x,2))
  for i in 1:size(x,2)
    jmax = first(remaining)
    for j in Iterators.drop(remaining,1)
      if matches[j,i] > matches[jmax,i]
        jmax = j
      end
    end

    delete!(remaining,jmax)
    order[i] = jmax
  end

  x[:,order]
end

function (tc::TCAnalysis)(x)
  frame_dims = size(x,2:ndims(x)...)
  frame_indices = collect(CartesianRange(frame_dims))[:]

  y = similar(x,div(size(x,1),tc.frame_len)+1,tc.N)
  N = prod(frame_dims)
  C = spzeros(eltype(x),N,N)

  Csparsity = 1-(1-tc.sparsity)^2
  c_t = Δt(tc.cort) / tc.τ

  old_ϕ = I
  @showprogress for t in 1:size(x,1)
    x_t = sparsify(@views(x[t,frame_indices]),tc.sparsity)
    C = C .+ c_t.*(x_t .* x_t' .- C) # note: .+= is slow (memory inefficient)

    if t % tc.frame_len == 0
      # find "objects" (principle components)
      _, ϕ = eigs(C,nev=tc.N)
      # maintain "object" order across time slices
      old_ϕ = ϕ = matchorder(ϕ,old_ϕ)

      # identify promenence of each object in the scene
      y[div(t,tc.frame_len)+1,:] .= abs.(ϕ'x_t)

      sparsify!(C,Csparsity)
    end
  end
  y
end

function rplot(tc::TCAnalysis,data::Matrix)
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)
  data = flipdim(data,2)

  df = DataFrame(response = data[:],
                 time = times(tc,data)[at(ixs,1)][:],
                 object_index = at(ixs,2)[:])

R"""

  library(ggplot2)

  ggplot($df,aes(x=time,y=object_index,fill=response)) +
    geom_raster() +
    ylab('Object Promenence') + xlab('Time (s)') +
    scale_fill_distiller(palette='Reds',name='Object ID',
                         direction=1)

"""
end
