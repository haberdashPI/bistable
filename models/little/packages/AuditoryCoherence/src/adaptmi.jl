using Parameters
using Unitful: s, ms, ustrip

# question: should we somehow limit the unit response???
# e.g. with:
sig(x) = 1/(1+exp(-10(x-0.5)))

@with_kw struct AdaptMI{S,I}
  α::Float64 = 1.0
  τ_y::typeof(1.0s) = 10ms
  shape_y::S = identity

  c_a::Float64 = 5
  τ_a::typeof(1.0s) = 1.5s

  c_e::Float64 = 0
  τ_e::typeof(1.0s) = 300ms

  c_d::Float64 = 0
  τ_d::typeof(1.0s) = 2s

  c_m::Float64 = 10
  τ_m::typeof(1.0s) = 50ms
  W_m::I = inhibit_uniform

  τ_σ::typeof(1.0s) = 100ms
  c_σ::Float64 = 0.3
end

# efficient representation of W*y where W = (ones(n,n) - I)/(n-1) for all n
function inhibit_uniform(y)
  x = (sum(y) .- vec(y)) ./ max(1,length(y)-1)
  reshape(x,size(y)...)
end

########################################
# the generic interface a type must implement for adaptmi to work on it

# empty_timeslice(x) = a single empty time slice of the object x, where x
# consists of multiple time slices
function empty_timeslice
end

# time_indices(x) = the number of time indices for x
function time_indices
end

# after set_timeslice!(y,t,y_t), time slice t of y === y_t
function set_timeslice!
end

# get_timeslice(y,t) is the t^th time slice of y
function get_timeslice
end

# __approx(f,x,xs...) ≈ f(x,xs...)

# x can be an 'important' value, meaning it is used to determine how the
# remaining values are approximated
function __approx
end

# get a timeslice in the same format used for values
# within the function called by __approx
function approx_empty_timeslice
end

# get a value with a similar format to data, but differing in that it has the
# same format as used within the function called by __approx
function approx_similar
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
approx_zeros(x,dims...) = zeros(x,dims...)

########################################
# complex number implementation
# - we 'approximate' values by their mangitude
#   so that m and a remain real, the resulting output
#   can be multiplied by the phase of the original signal
#   to find the final, complex output
function __approx(f,x::AbstractArray{<:Complex},args...)
  y = f(abs.(x),args...)
  withangle(angle.(x),y)
end
function withangle(angle,y::Tuple)
  # @show angle
  # @show y[1]
  (y[1].*exp.(angle.*im),y[2:end]...)
end
function withangle(angle,y)
  y.*exp.(angle.*im)
end
approx_similar(x::AbstractArray{<:Complex{T}}) where T =
  similar(x,T)
approx_zeros(y::AbstractArray{<:Complex{T}}) where T =
  zeros(eltype(y),size(y)...)
approx_zeros(y::AbstractArray{<:Complex{T}},dims...) where T =
  zeros(eltype(y),dims...)

################################################################################
# genertic adaptation and mutual-inhibition operation
adaptmi(x;progressbar=true,kw...) = adaptmi(x,AdaptMI(;kw...),progressbar)
function adaptmi(x,params::AdaptMI,progressbar=true)
  @assert :time ∈ axisnames(x)
  time = Axis{:time}

  y = similar(x)

  α = params.α
  τ_y = params.τ_y; shape_y = params.shape_y
  τ_a = params.τ_a; c_a = params.c_a
  τ_e = params.τ_e; c_e = params.c_e
  τ_d = params.τ_d; c_d = params.c_d
  τ_m = params.τ_m; c_m = params.c_m; W_m = params.W_m

  a = approx_similar(y) # a = adaptation
  e = approx_similar(y) # a = adaptation
  d = approx_similar(y) # a = adaptation
  m = approx_similar(y) # m = mutual inhibition

  yr_t = zeros(y[time(1)])
  a_t = approx_zeros(yr_t)
  e_t = approx_zeros(yr_t)
  d_t = approx_zeros(yr_t)
  m_t = approx_zeros(yr_t)

  dt_y = eltype(a_t)(Δt(x) / τ_y)
  dt_a = eltype(a_t)(Δt(x) / τ_a)
  dt_e = eltype(e_t)(Δt(x) / τ_e)
  dt_d = eltype(d_t)(Δt(x) / τ_d)
  dt_m = eltype(a_t)(Δt(x) / τ_m)

  progress = progressbar ? Progress(desc="Adapt/Inhibit: ",ntimes(y)) : nothing
  for t in indices(times(y),1)
    y_t = x[time(t)]

    yr_t,a_t,e_t,d_t,m_t = @approx (y_t,a_t,e_t,d_t,m_t) begin
      y_t .*= α
      y_t .-= (y_t.*c_a.*a_t .- c_e.*e_t .+ c_d.*d_t .+ c_m.*m_t)
      yp_t = shape_y.(y_t)
      a_t .+= (yp_t .- a_t).*dt_a
      e_t .+= (yp_t .- e_t).*dt_e
      d_t .+= (yp_t .- d_t).*dt_d
      m_t .+= (W_m(yp_t) .- m_t).*dt_m

      yp_t,a_t,e_t,d_t,m_t
    end

    y[time(t)] = yr_t
    a[time(t)] = a_t
    e[time(t)] = e_t
    d[time(t)] = d_t
    m[time(t)] = m_t

    next!(progress)
  end
  y,a,m,e,d
end

################################################################################
# generic drifting noise function

drift(x,along_axes...;progressbar=true,kw...) =
  drift(x,AdaptMI(;kw...),along_axes,progressbar)
function drift(x,params::AdaptMI,along_axes=typeof.(axes(x)),progressbar=true)
  τ_σ, c_σ = params.τ_σ, params.c_σ
  time = Axis{:time}

  y = similar(x)
  temp = zeros(x[time(1),(ax(1) for ax in along_axes)...])
  σ_t = approx_zeros(temp)
  dims = size(σ_t)
  progress = progressbar ? Progress(desc="Drift: ",ntimes(y)) : nothing
  for t in indices(times(y),1)
    y_t = x[time(t)]
    y_t,σ_t = @approx (y_t,σ_t) begin
      σ_t .+= -σ_t.*(Δt(x)/τ_σ) .+ randn(dims).*(c_σ*sqrt(2Δt(x)/τ_σ))
      y_t .*= (1 .+ σ_t)

      y_t,σ_t
    end

    y[time(t)] = y_t

    next!(progress)
  end
  y
end
