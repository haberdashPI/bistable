# TODO: for updating via adaptation I think it should be possible to think in
# terms of the updates happening to the PCA representation, and then show the
# equivalence of that to the updates happening to the covariance matrix
# (and/or, do we care if it is equivalent, perhaps this is the representation
# we want to operate over)

using Parameters
include("cortical.jl")
include("online_pca.jl")

abstract type TCAnalysis end
Δt(tc::TCAnalysis) = Δt(tc.upstream)
times(tc::TCAnalysis,x) = times(tc.upstream,x)

struct OnlineTCAnalysis <: TCAnalysis
  upstream::CorticalModel
  ncomponents::Int
  method::Symbol
  init_len::Int
end

(tc::OnlineTCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))

# alternative: I could have a different
# set of eigenvectors for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (tc::OnlineTCAnalysis)(x)
  tdims = size(x,2)
  fdims = size(x,3,4)

  tI = collect(CartesianRange(tdims))
  fI = collect(CartesianRange(fdims))

  λ = similar(x,size(x,1),tc.ncomponents)
  y = similar(x,size(x,1),tc.ncomponents)
  ϕ = similar(x,size(x,1),tc.ncomponents,fdims...)
  λ[1:start-1,:] = 0
  ϕ[1:start-1,:] = 0
  y[1:start-1,:] = 0

  for t in start:size(x,1)
    x_t = reshape(view(x,t,:,:,:),size(x,2),:)
    for ti in size(x_t,1)
      C = LinAlg.lowrankupdate(C,x_t[ti,:])
    end
    λ[t,:] = eigvals(C)
    ϕ[t,:,fI] = reshape(eigvecs(C)',:,fdims...)
    y[t,:] = mean(x_t*eigvecs(C),1)
  end

  λ,ϕ,y
end

struct BatchTCAnalysis <: TCAnalysis
  upstream::CorticalModel
  ncomponents::Int
end

(tc::BatchTCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
function (tc::BatchTCAnalysis)(x)
  sv, = svds(reshape(x,prod(size(x,1,2)),:),nsv=tc.ncomponents)

  yv = reshape(x,prod(size(x,1,2)),:)*sv[:V]
  y = reshape(yv,size(x,1:2...)...,:)
  return sv[:S].^2 / size(x,1), reshape(sv[:V]',:,size(x,3,4)...), y
end

TCAnalysis(upstream,ncomponents;method=:ipca,init_len=0) =
  method ∈ [:batch,:pca] ? BatchTCAnalysis(upstream,ncomponents) :
  OnlineTCAnalysis(upstream,ncomponents,method,init_len)


#========================================
struct TCAnalysis
  upstream::CorticalModel
  ncomponents::Int
end
times(tc::TCAnalysis,x) = times(tc.upstream,x)
times(tc::TCAnalysis,x::Array{T,3}) where T = times(tc.upstream,x[:,1,:])
freqs(tc::TCAnalysis,x::Array{T,3}) where T = freqs(tc.upstream,x[:,1,:])

(tc::TCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
function (tc::TCAnalysis)(x)
  sv, = svds(reshape(x,prod(size(x,1,2)),:),nsv=tc.ncomponents)

  sv[:S].^2 / size(x,1), reshape(sv[:V],:,size(x,3,4))
end
========================================#

function rplot(tc::TCAnalysis,data::Matrix)
  ixs = CartesianRange(size(abs.(data)))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(times(tc,data)[at(ixs,1)]),
                 component_index = vec(at(ixs,2)))

R"""

  library(ggplot2)

  ggplot($df,aes(x=time,y=component_index,fill=response)) +
    geom_raster() +
    ylab('Component') + ('Time (s)') +
    scale_fill_distiller(palette='Reds',name='Amplitude',
                         direction=1)

"""
end

function rplot(tc::TCAnalysis,data::Array{T,3}) where T
  ixs = CartesianRange(size(abs.(data)))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(times(tc,data)[at(ixs,1)]),
                 component_index = vec(at(ixs,2)),
                 freq_index = vec(at(ixs,3)))

  fbreaks = 2.0.^(-3:2)
  fs = freqs(tc,data)
  findices = mapslices(abs.(1000.0.*fbreaks .- fs'),2) do row
    _, i = findmin(row)
    i
  end

R"""

  library(ggplot2)

  ggplot($df,aes(x=time,y=component_index,fill=response)) +
    geom_raster() +
    facet_wrap(~component_index) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz') + xlab('Time (s)') +
    scale_fill_distiller(palette='Reds',name='Amplitude',
                         direction=1)

"""
end
