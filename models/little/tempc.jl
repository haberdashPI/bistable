# TODO: for updating via adaptation I think it should be possible to think in
# terms of the updates happening to the PCA representation, and then show the
# equivalence of that to the updates happening to the covariance matrix
# (and/or, do we care if it is equivalent, perhaps this is the representation
# we want to operate over)

using Parameters
include("cortical.jl")
include("online_pca.jl")
#=========================================
# TODO: allow TC anlaysis to usee batch PCA
# method

abstract type TCAnalysis end
Δt(tc::TCAnalysis) = Δt(tc.upstream)
times(tc::TCAnalysis,x) = times(tc.upstream,x)

struct OnlineTCAnalysis <: TCAnalysis
  upstream
  ncomponents::Int
  method::Symbol
  init_len::Int
end

(tc::OnlineTCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))

Base.CartesianRange(x::Int) = CartesianRange((x,))
function (tc::OnlineTCAnalysis)(x)
  ii = vec(collect(CartesianRange(size(x,2:ndims(x)...))))
  if tc.init_len > 0
    sv, = svds(x[1:tc.init_len,ii],nsv=tc.ncomponents)
    C = OnlinePCA(sv,tc.init_len,method=tc.method)
    start = tc.init_len+1
  else
    C = OnlinePCA(length(ii),tc.ncomponents,method=tc.method)
    start = 1
  end

  y = similar(x,size(x,1),tc.ncomponents)
  y[1:start-1,:] = 0
  for t in start:size(x,1)
    x_t = x[t,ii]
    C = LinAlg.lowrankupdate(C,x_t)
    y[t,:] = eigvecs(C)' * x_t
  end

  y[:,sortperm(eigvals(C))],C
end

struct BatchTCAnalysis <: TCAnalysis
  upstream
  ncomponents::Int
end

(tc::BatchTCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
function (tc::BatchTCAnalysis)(x)
  ii = vec(collect(CartesianRange(size(x,2:ndims(x)...))))
  @views begin
    sv, = svds(x[:,ii],nsv=tc.ncomponents)

    x[:,ii]*sv[:V][:,sortperm(sv[:S])],sv
  end
end

TCAnalysis(upstream,ncomponents;method=:ipca,init_len=0) =
  method ∈ [:batch,:pca] ? BatchTCAnalysis(upstream,ncomponents) :
  OnlineTCAnalysis(upstream,ncomponents,method,init_len)
=========================================#

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

  sv[:S].^2 / size(x,1), sv[:V]
end

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
    ylab('Object ID') + xlab('Time (s)') +
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
