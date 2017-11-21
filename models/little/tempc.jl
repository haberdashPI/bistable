using Parameters
include("cortical.jl")

struct SparseCorrelationRange
  I::Vector{Int}
end
function length(x::SparseCorrelationRange)
  n = length(x.I)
  n + (n-1)n>>1
end

start(x::SparseCorrelationRange) = (1,1)

function next(x::SparseCorrelationRange,itr::NTuple{2,Int})
  i,j = itr
  result = CartesianIndex(i,j)

  n = length(x.I)
  if j < n
    result,(i,j+1)
  else
    result,(i+1,i+1)
  end
end

function done(x::SparseCorrelationRange,itr::NTuple{2,Int})
  i,j = itr
  n = length(x.I)

  i >= n && j >= n
end

mutable struct SparseCorrelation
  vals::Vector{Int}
  II::Vector{CartesianIndex{2}}
  len::Int
  N::Int
end

function SparseCorrelation(N,sparsity)
  # make a guess about how much space we'll need
  n = ((1-sparsity)N)^1.5
  SparseCorrelation(Array{Int}(n),Array{CartesianIndex{2}}(n),0,N)
end

function update!(C::SparseCorrelation,n)
  if n > length(C.vals)
    resize!(C.vals,n)
    resize!(C.II,n)
  end
  C.len = n
end

indices(C::SparseCorrelation) = filter(ii -> ii[1] > 0,C.II)

function sorted_merge(X,Y,lt,eq,XT,YT)
  result = CircularDeque{Union{XT,YT}}(length(X) + length(Y))
  state1,state2 = start(X),start(Y)
  item1 = zero(XT)
  item2 = zero(YT)

  sorted_merge_helper!(result,X,Y,lt,eq,state1,state2,item1,item2)
end

function sorted_merge_helper!(result,X,Y,lt,eq,state1,state2,item1,item2)
  a1,a2 = false,false
  d1,d2 = done(state1),done(state2)
  while !((d1 && !a1) && (d2 && !a2))
    if !d1 && !a1 && !d2 && !a2
      item1,state1 = next(X,state1)
      item2,state2 = next(X,state2)
      a1 = a2 = true
      d1,d2 = done(state1),done(state2)
    elseif !d2 && !a2
      item2,state2 = next(X,state2)
      a2 = true
      d2 = done(state2)
    elseif !d1 && !a1
      item1,state1 = next(X,state2)
      a1 = true
      d1 = done(state1)
    end

    if lt(item1,item2)
      a1 = false
      push!(result,item1)
    elseif eq(item1,item2)
      a1 = a2 = false
      push!(result,item1)
    else
      a2 = false
      push!(result,item2)
    end
  end

  result
end

function update_correlation(C::SparseCorrelation,x,θ,c_t)
  large = find(x -> abs(x) > θ,x)
  xabs = [abs(x[i]) for i in large]
  maxabs = maximum(xabs)
  entries = nonzero[find(i -> maxabs*xabs[i] > C.θ,nonzero)]

  newi = SparseCorrelationRange(entries)
  oldi = filter(ii -> ii[1] > 0,1:C.len)

  indices = sorted_merge(newi,oldi,(i,j) -> corless(newi,C.II[j]),
                         (i,j) -> i == C.II[j])

  update!(w,length(indices))

  @threads for (k,II) in enumerate(indices)
    i,j = II.T
    v = c_t*(x[i] * x[j] - C[i,j])

    if v > C.θ
      w.vals[k] = v
      w.I[k] = i
      w.J[k] = j
    else
      w.I[k] = 0
    end
  end
end

@with_kw struct TCAnalysis
  cort::CorticalModel
  τ::Seconds{Float64} = 1s
  sparsity::Float64 = 0.8
end
Δt(tc::TCAnalysis) = Δt(tc.cort)

# effecient algorithm to find the quantile of a sparse matrix
# when we know there are no negative entries
function aquantile(x::AbstractSparseMatrix,p)
  @assert 0 <= p <= 1 "p must be in the interval [0,1]"
  _,_,nz = findnz(x)
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

# NOTE: there is potentially a real problem with sparsity and gradual changes in
# C.
function (tc::TCAnalysis)(x;drop_interval=10)
  y = similar(x)
  N = prod(size(x)[[1,2,4]])
  C = spzeros(eltype(x),N,N)
  c_t = Δt(tc) / tc.τ

  φ = 1-(1-tc.sparsity)^2
  drop_wait = 0
  C = SparseCorrelation()
  for t in 1:size(x,3)
    update_correlation(C,x,φ,c_t)
    @views y[:,:,t,:][:] = C*x[:,:,t,:][:]

    aC = abs.(C)
    cut = aquantile(aC,φ)
    iis = find(0 .< aC .< cut)
    C[iis] = 0
    if length(iis) > 0
      if (drop_wait += 1) > drop_interval
        drop_wait = 0
        dropzeros!(C)
      end
    end
  end
  y
end



# struct CorMerge{T}
#   merged::T
# end
# merge_correlation_indices(x,y) = CorMerge((x,y))
# function corless(x::CartesianIndex{2},y::CartesianIndex{2})
#   x[1] < y[1] ? true :
#   x[1] > y[1] ? false :
#   x[2] < y[2]
# end

# struct CorMergeItr{T}
#   status::Int
#   item::NTuple{2,CartesianIndex{2}}
#   data::T
# end

# function start(x::CorMerge)
#   state1,state2 = start(x.merged[1]),start(x.merged[2])
#   CorMergeItr(-done(x.merged[1],state1) - 2done(x.merged[2],state2),
#               (CartesianIndex(0,0),CartesianIndex(0,0)),
#               (state1,state2))
# end

# function next(x::CorMerge,state)
#   function helper(item1,item2,state1,state2)
#     if corless(item1,item2)
#       if done(x.merged[1],state1)
#         item1,CorMergeItr(-1,(item1,item2),(state1,state2))
#       else
#         item1,CorMergeItr(2,(item1,item2),(state1,state2))
#       end
#     elseif state[1] == state[2]
#       item1,CorMergeItr(-done(x.merged[1],state1) - 2done(x.merged[2],state2),
#                         (item1,item2),(state1,state2))
#     else
#       item2,CorMergeItr(1,(item1,item2),(state1,state2))
#     end
#   end

#   if state.status == 0
#     item1,state1 = next(x.merged[1],state.data[1])
#     item2,state2 = next(x.merged[2],state.data[2])
#     helper(item1,state1,item2,state2)
#   elseif state.status == 2
#     item1,state1 = next(x.merged[1],state.data[1])
#     item2,state2 = state.item[2],state.data[2]
#     helper(item1,state1,item2,state2)
#   elseif state.status == 1
#     item1,state1 = state.item[1],state.data[1]
#     item2,state2 = next(x.merged[2],state.data[2])
#     helper(item1,state1,item2,state2)
#   elseif state.status == -1
#     item2,state2 = next(x.merged[2],state.data[2])
#     item1,state1 = state.item[1],state.data[1]
#     item2,CorMergeItr(-1 - 2done(x.merged[2],item2),
#                       (item1,item2),(state1,state2))
#   elseif state.status == -2
#     item1,state1 = next(x.merged[1],state.data[1])
#     item2,state2 = state.item[2],state.data[2]
#     item1,CorMergeItr(-2 -done(x.merged[1],item1),
#                       (item1,item2),(state1,state2))
#   else
#     error("No items available")
#   end
# end
# done(x::CorMerge) = x.status == -3

# this is an optimized version of the following operations
#
# C .+= c_t.*(x .* x' .- C)
# C[abs.(C) .<= θ] = 0
