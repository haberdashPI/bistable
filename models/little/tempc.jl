using Parameters
using ProgressMeter
include("cortical.jl")

@with_kw struct TCAnalysis <: CorticalRepresentation
  cort::CorticalModel
  τ::Seconds{Float64} = 1s
  sparsity::Float64 = 0.8
  frame_len = 10
  prior::Float64 = 1.0
end
Δt(tc::TCAnalysis) = Δt(tc.cort) * tc.frame_len

scales(tc::TCAnalysis) = scales(tc.cort)
rates(tc::TCAnalysis) = rates(tc.cort)
freqs(tc::TCAnalysis,y) = freqs(tc.cort,y)
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
end

function normalize(x_i,covar,shift=1.0)
  # multidimensional equivalent of x_i / sd(x)

  # NOTE: the somewhat akward spelling here
  # is due to the limited number of methods implemented
  # for sparse matrices
  cholfact(covar,shift=shift)[:UP] \ x_i[:,:]
end

function (tc::TCAnalysis)(x)
  frame_dims = size(x,2:4...)
  frame_indices = collect(CartesianRange(frame_dims))[:]

  y = zeros(typeof(x[1]),div(size(x,1),tc.frame_len)+1,frame_dims...)
  N = prod(frame_dims)
  C = spzeros(eltype(x),N,N)
  Csparsity = 1-(1-tc.sparsity)^2
  c_t = Δt(tc.cort) / tc.τ

  sparsity = 0
  total = 0
  @showprogress for t in 1:size(x,1)
    x_t = sparsify(@views(x[t,frame_indices]),tc.sparsity)
    C = C .+ c_t.*(x_t .* x_t' .- C) # note: .+= is slow (memory inefficient)

    if t % tc.frame_len == 0
      C = sparsify!(C,Csparsity)
      xnorm = normalize(x_t,C,tc.prior)
      y[div(t,tc.frame_len)+1,frame_indices] = xnorm

      sparsity += length(xnorm.nzval)
      total += length(xnorm)
    end
  end

  per = round(100(1-(sparsity / total)),1)
  perstr = per > 99.9 ? ">99.9" : string(per)
  println("Achieved an average sparsity of $(per)%")
  y
end
