using Combinatorics

mutable struct MultiNormalStats{T}
  x::Vector{T}
  x2::Matrix{T}
  n::Int
  x2_offset::Int
end
function MultiNormalStats(mean,cov,n::Int,x2_offset=length(mean))
  MultiNormalStats(mean*n,cov .+ mean.*mean',n,x2_offset)
end
function MultiNormalStats(data::AbstractArray{T,3} where T,
                          n=size(data,1),x2_offset=size(data,2))
  x = similar(data,size(data,2))
  x2 = similar(data,size(data,2,2))
  for k in size(data,3)
    d = data[:,:,k]
    x .+= vec(mean(d,1))*n
    x2 .+= n .* d'd ./ size(data,1)
  end
  MultiNormalStats(x,x2,n,x2_offset)
end

@with_kw struct PriorTracking <: Tracking
  tc::typeof(1.0s) = 1s
  prior::MultiNormalStats{Float64}
end
Tracking(::Val{:prior};params...) = PriorTracking(;params...)

function track(C::Coherence,params::PriorTracking)
