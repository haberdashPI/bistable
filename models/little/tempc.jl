using Parameters
using ProgressMeter
include("cortical.jl")

@with_kw struct TCAnalysis
  cort::CorticalModel
  τ::Seconds{Float64} = 1s
  sparsity::Float64 = 0.8
end
Δt(tc::TCAnalysis) = Δt(tc.cort)

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

function sparsify(x::AbstractVector,θ)
  xa = abs.(x)
  cut = quantile(xa,θ)
  indices = find(x -> x > cut,xa)
  SparseVector(length(x),indices,x[indices])
end

const num_dropped = fill(0)
const drop_min = 4 * 2^10
function sparsify(C::SparseMatrixCSC,θ)
  aC = abs.(C)
  cut = aquantile(aC,θ)
  iis = find(0 .< aC.nzval .< cut)
  if length(iis) > 0
    C.nzval[iis] = 0
    if (num_dropped[] += length(iis)) > drop_min
      num_dropped[] = 0
      dropzeros!(C)
    end
  end
end

function (tc::TCAnalysis)(x)
  y = similar(x)
  N = prod(size(x)[[1,2,4]])
  C = spzeros(eltype(x),N,N)
  c_t = Δt(tc) / tc.τ

  Csparsity = 1-(1-tc.sparsity)^2
  C = spzeros(N,N)
  sparsity = 0
  @showprogress for t in 1:size(x,3)
    xi = sparsify(@views(x[:,:,t,:][:]),tc.sparsity)
    # NOTE: we cannot use .+= becuase this leads to problematic memory usage
    C = C .+ c_t.*(xi .* xi' .- C)
    xn = C*xi
    @views y[:,:,t,:][:] = xn

    sparsify(C,Csparsity)
    sparsity += length(xn.nzval)
  end
  per = round(100(1-(sparsity / prod(size(x)))))
  perstr = per > 99 ? ">99" : string(per)
  println("Achieved an average sparsity of $(per)%")
  y
end
