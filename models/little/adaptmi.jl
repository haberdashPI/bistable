using Parameters
using Unitful: s, ms, ustrip

include("online_pca.jl")
# question: should we somehow limit the unit response???
# e.g. with:
sig(x) = 1/(1+exp(-10(x-0.5)))

# TODO: because of some of the performance issues with keyword arguments, I use
# a params object.  However, once julia v1.0 hits, I should be able to get
# similar performance using just keyword arguments
@with_kw struct AdaptMI{S,I}
  τ_y::Seconds{Float64} = 10ms
  shape_y::S = identity

  c_a::Float64 = 5
  τ_a::Seconds{Float64} = 1.5s

  c_m::Float64 = 10
  τ_m::Seconds{Float64} = 50ms
  W_m::I = inhibit_uniform

  Δt::Seconds{Float64} = 1ms
end

# efficient representation of W*y where W = (ones(n,n) - I)/(n-1) for all n
function inhibit_uniform(y)
  x = (sum(y) .- vec(y)) ./ max(1,length(y)-1)
  reshape(x,size(y)...)
end

# by default a call to adaptmi updates y_t at each time step so that it
# incrementally approaches the current value of x_t
function adaptmi(x::AbstractArray{T},params=AdaptMI()) where T
  adaptmi(similar(x),params) do y_t,t,dt
    I = CartesianRange(size(y_t))
    for ii in I; y_t[ii] += (x[t,ii] - y_t[ii])*dt; end
    y_t
  end
end

function noisey!(x,τ_ε,c_ε,Δt=1ms)
  ε_t = zeros(size(x,2))
  n = size(x,2)
  for t in 1:size(x,1)
    ε_t .+= -ε_t.*(Δt/τ_ε) + randn(n).*(c_ε*sqrt(2Δt/τ_ε))
    x[t,:] += ε_t
  end
  x
end
noisey(x,τ_ε,c_ε) = noisey!(copy(x),τ_ε,c_ε)


########################################
# the generic interface a type must implement for adaptmi to work on it

# timeslice(x) = a single empty time slice of the object x, where x
# consists of multiple time slices
function timeslice
end

# nunits(x_t) =  the number of units (dimensions) for a single time slice
function nunits
end

# time_indices(x) = the number of time indices for x
function time_indices
end

# after set_timeslice!(y,t,y_t), time slice t of y === y_t
function set_timeslice!
end

# __approx(f,x,xs...) ≈ f(x,xs...)

# x can be an 'important' value, meaning it is used to determine how the
# remaining values are approximated
function __approx
end

# @approx (x1,x2) begin
#   x1 + x2
# end
#
# is equivalent to:
#
# __approx(+,x1,x2)
macro approx(args,body)
  :(__approx($(esc(args)) -> $(esc(body)),$(esc(args))...))
end

##############################
# default implementation
# - we assume an array like object and treat the first dimension
#   as the time dimension
# - there is no loss due to approximation

__approx(f,args...) = f(args...)
approx_similar(x) = similar(x)
timeslice(y) = zeros(eltype(y),size(y)[2:end]...)
approx_timeslice(y) = timeslice(y)
nunits(y_t) = length(y_t)
time_indices(y) = indices(y,1)
function set_timeslice!(y,t,y_t)
  for ii in CartesianRange(size(y_t)); y[t,ii] = y_t[ii]; end
  y_t
end

########################################
# complex number implementation
# - we 'approximate' values by their mangitude
#   so that m and a remain real, the resulting output
#   can be multiplied by the phase of the original signal
#   to find the final, complex output
function __approx(f,x::Array{<:Complex},args...)
  y = f(abs.(x),args...)
  withangle(angle.(x),y)
end
function withangle(angle,y::Tuple)
  (y[1].*angle,y[2:end]...)
end
withangle(angle,y) = y.*angle.(x)
approx_similar(x::Array{<:Complex{T}}) where T = similar(x,T)
approx_timeslice(y::Array{<:Complex{T}}) where T = zeros(T,size(y)[2:end]...)

##############################
# implementation for EigenSeries
# - each time slice is an EigenSpace
# - to approximate f(x,ys...)  we project all secondary eigenspaces ys onto the
#   eigenspace x and then perform the operation over the resulting eigenvalues
timeslice(y::EigenSeries) = EigenSpace(y)
nunits(x::EigenSpace) = ncomponents(x)
time_indices(x::EigenSeries) = indices(x.u,1)
set_timeslice!(x::EigenSeries,t,slice) = x[t] = slice

function __approx(f,x::EigenSpace,ys::EigenSpace...)
  λ = f(x.λ,(project(x,y).λ for y in ys)...)
  eigenspaces(x,λ)
end
# transform eigenvalues into eigenspace
eigenspaces(x,λ::Vector) = EigenSpace(x.u,λ)
# transform tuple of eigenvalues into tuple of eigenspaces
eigenspaces(x,λs::Tuple) = map(λ -> EigenSpace(x.u,λ),λs)

################################################################################
# genertic adaptation and mutual-inhibition operation
function adaptmi(update,y,params)
  τ_y = params.τ_y; shape_y = params.shape_y
  τ_a = params.τ_a; c_a = params.c_a
  τ_m = params.τ_m; c_m = params.c_m; W_m = params.W_m
  Δt = params.Δt

  a = approx_similar(y) # a = adaptation
  m = approx_similar(y) # m = mutual inhibition

  yr_t = timeslice(y)
  a_t = approx_timeslice(y)
  m_t = approx_timeslice(y)

  dt_y = Δt / τ_y
  dt_a = Δt / τ_a
  dt_m = Δt / τ_m

  @showprogress "Adapt/Inhibit: " for t in time_indices(y)
    y_t = update(yr_t,t,dt_y)

    yr_t,a_t,m_t = @approx (y_t,a_t,m_t) begin
      y_t .-= (y_t.*c_a.*a_t .+ c_m.*m_t)
      yp_t = shape_y.(y_t)
      a_t .+= (yp_t .- a_t).*dt_a
      m_t .+= (W_m(yp_t) .- m_t).*dt_m

      y_t,a_t,m_t
    end

    set_timeslice!(y,t,yr_t)
    set_timeslice!(a,t,a_t)
    set_timeslice!(m,t,m_t)
  end
  y,a,m
end
