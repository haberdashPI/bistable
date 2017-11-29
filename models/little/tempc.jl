using Parameters
using ProgressMeter
include("cortical.jl")
include("online_pca.jl")

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

  y[:,sortperm(eigvals(C))]
end

struct BatchTCAnalysis <: TCAnalysis
  upstream
  ncomponents::Int
end

function (tc::BatchTCAnalysis)(x)
  ii = vec(collect(CartesianRange(size(x,2:ndims(x)...))))
  @views begin
    sv, = svds(x[:,ii],nsv=tc.ncomponents)

    x[:,ii]*sv[:V][:,sortperm(sv[:S])]
  end
end

TCAnalysis(upstream,ncomponents;method=:ipca,init_len=0) =
  method ∈ [:batch,:pca] ? BatchTCAnalysis(upstream,ncomponents) :
  OnlineTCAnalysis(upstream,ncomponents,method,init_len)

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
    scale_fill_distiller(palette='Reds',name='Object "Strength"',
                         direction=1)

"""
end
