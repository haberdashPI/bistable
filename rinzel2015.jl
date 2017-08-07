using Plots
using Parameters
using Base.Iterators: drop, take, partition
using Lazy

@with_kw type BistableParams
  # threshold parameters
  θ_F::Float64 = 0.2  # threshold
  k_F::Float64 = 12   # threshold slope

  # input parameters
  Λ::Float64 = 1/6      # plateau fraction
  α_1::Float64 = 0.015  # impulse rise time
  α_2::Float64 = 0.0825 # plateau stay time
  I_p::Float64 = 0.525  # amplitude
  σ_p::Float64 = 8      # input decy in semitones

  # intrinsic parametesr
  g::Float64 = 0.065  # adaptation strength
  γ::Float64 = 0.075  # noise strength
  β_i::Float64 = 0.3  # inhibition strength
  σ_i::Float64 = 10   # inhibition decay constant
  β_e::Float64 = 0.7  # recurrent excitation
  κ::Float64 = 0      # excitation decay strength

  # timescales
  τ_r::Float64 = 0.010 # 'cortical' timescale (change in rate of firing)
  τ_a::Float64 = 1.4   # adaptation timescale
  τ_e::Float64 = 0.070 # NMDA excitation timescale
  τ_χ::Float64 = 0.100 # noise timescale
  τ_d::Float64 = NaN   # not present in first model
end

const second_model_params = BistableParams(
  I_p = 0.47,
  σ_p = 8.5,
  σ_i = Inf,
  β_e = 0.85,
  κ = 0.25,
  τ_d = 3.0
)

const n_params = 18
const n_noise = 3

const r = [1,2,3]
const a = [4,5,6]
const e = [7,8,9]
const d = [10,11,12]
const χ = [13,14,15]
const Σ_I = [16,17,18]

type Simulation
  t::Range{Float64}
  us::Base.Generator
  If::AbstractArray{Float64,2}
end

function OfflineInputs(DF,PR,sim_length;dt=1/1000,params=BistableParams())
  const p = params
  const t = 0:dt:sim_length
  const f_Is = Float64[1,1+DF]
  const f = Float64[1,1+DF/2,1+DF]

  triplets = floor(Int,PR*sim_length/4)
  aba_fs = repeat([1,2,1],outer=triplets)
  aba_ts = cumsum([0; repeat([1/PR,1/PR,2/PR],outer=triplets)])[1:end-1]

  function input(t::Float64)
    if t > 0.0
      (exp(2)/p.α_1^2) * t^2*exp(-2t / p.α_1) +
        (p.Λ*exp(2)/p.α_2^2) * t^2*exp(-2t / p.α_2)
    else 0.0
    end
  end

  If = zeros(length(t),2)
  max_t_len = floor(Int,PR*10/dt)
  for i in eachindex(aba_fs)
    fi = aba_fs[i] == f_Is[1] ? 1 : 2

    start_t_i = floor(Int,aba_ts[i]/dt) + 1
    for t_i in start_t_i:(start_t_i + max_t_len)
      t_i <= length(t) || break
      If[t_i,fi] += input(t[t_i] - aba_ts[i])
    end
  end

  If
end

type OnlineInputs <: AbstractArray{Float64,2}
  freqis::Vector{Float64}
  times::Vector{Float64}
  t::Range{Float64}
  f::Function
  max_t::Float64
end

function OnlineInputs(DF,PR,sim_length;dt=1/1000,params=BistableParams())
  const p = params
  t = 0:dt:sim_length
  triplets = floor(Int,PR*sim_length/4)

  function input(t::Float64)
    if t > 0.0
      (exp(2)/p.α_1^2) * t^2*exp(-2t / p.α_1) +
        (p.Λ*exp(2)/p.α_2^2) * t^2*exp(-2t / p.α_2)
    else 0.0
    end
  end

  freqis = repeat([1,2,1],outer=triplets)
  times = cumsum([0; repeat([1/PR,1/PR,2/PR],outer=triplets)])[1:end-1]
  max_t = 10*1/PR

  OnlineInputs(freqis,times,t,input,max_t)
end

Base.size(x::OnlineInputs) = (length(x.t),2)
function Base.getindex(x::OnlineInputs,i::Int,j::Int)
  ks = find((x.times .<= x.t[i] .< (x.times .+ x.max_t)) .& (x.freqis .== j))
  if isempty(ks) 0.0
  else sum(k -> x.f(x.t[i] - x.times[k]),ks) end
end

function Simulation(DF,PR,sim_length;dt=1/1000,
                    params=BistableParams(),inputs=OfflineInputs)
  If = inputs(DF,PR,sim_length,dt=dt,params=params)

  const p = params
  const t = 0:dt:sim_length
  const f_Is = Float64[1,1+DF]
  const f = Float64[1,1+DF/2,1+DF]

  # NOTE: the example inputs shown in the paper do not match the parameters for
  # w(Δf) as specificied in the paper.

  w(Δf::Float64) = p.I_p * exp(-Δf / p.σ_p)
  C_i(Δf::Float64) = p.β_i * exp(-Δf^2 / (2*p.σ_i^2))
  F(x::Float64) = 1 / (1 + exp(p.k_F*(-x + p.θ_F)))

  function bistable(i::Int,t::Float64,u_0::Vector{Float64})
    u_1 = copy(u_0)
    for k in 1:3
      Σ_r = sum(h -> C_i(abs(f[k] - f[h]))*u_0[r[h]],1:3)
      u_1[Σ_I[k]] = sum(j -> w(abs(f_Is[j] - f[k])) * If[i,j],1:2)

      u_1[r[k]] +=
        (-u_0[r[k]] + F(p.β_e*u_0[e[k]]*u_0[d[k]] - Σ_r - p.g*u_0[a[k]] +
                        u_1[Σ_I[k]] + u_0[χ[k]]))/p.τ_r * dt

      u_1[a[k]] += (-u_0[a[k]] + u_0[r[k]])/p.τ_a * dt
      u_1[e[k]] += (-u_0[e[k]] + u_0[r[k]])/p.τ_e * dt

      if isnan(p.τ_d)
        u_1[d[k]] = 1
      else
        u_1[d[k]] += (-u_0[d[k]] + (1-p.κ*u_0[r[k]]))/p.τ_d * dt
      end

      u_1[χ[k]] += -u_0[χ[k]] / p.τ_χ * dt +
        p.γ * √(2/p.τ_χ) * sqrt(dt) * randn()
    end

    u_1
  end

  u = 1e-6*rand(n_params)
  us = (u = bistable(i,t_i,u) for (i,t_i) in enumerate(t))

  Simulation(t,us,If)
end

function asmatrix(itr,N)
  x = first(itr)
  xs = Array{Float64,2}(N,length(x))
  xs[1,:] = x
  i = 1
  for x in drop(itr,1)
    xs[i+=1,:] = x
    i < size(xs,1) || break
  end

  xs
end

const τ_response = 0.05
function response_fn(sim)
  let dt = sim.t[2] - sim.t[1],
    α = dt / (τ_response + dt),
    y = 0

    function response(u)
      y = α*u[r] + (1-α)y
      y[2] > (y[1] + y[3])/2
    end
  end
end

function Base.collect(sim::Simulation)
  U = asmatrix(sim.us,length(sim.t))
  return sim.t,U,sim.If,mapslices(response_fn(sim),U,2)
end

function partial_collect(sim::Simulation;
                         min=minimum(sim.t),max=maximum(sim.t),decim=1)
  dt = sim.t[2] - sim.t[1]
  min_i = floor(Int,min / dt) + 1
  max_i = floor(Int,max / dt) + 1
  ii = min_i:decim:max_i

  partial = @> sim.us take(max_i) drop(min_i-1) partition(decim)
  U = asmatrix((first(xs) for xs in partial),div(max_i - min_i,decim))
  return sim.t[ii],U,sim.If[ii,:],mapslices(response_fn(sim),U,2)
end

function responses(sim::Simulation)
  f = response_fn(sim)
  return (f(u) for u in sim.us)
end


# OLD CODE:
#   rates = u[:,r]

#   rates[1,:] *= α
#   for i in 1:(size(rates,1)-1)
#     rates[i+1,:] = α * rates[i+1,:] + (1-α) * rates[i,:]
#   end
#   rates
# end


# function responses(t,u)
#   rs = filter_response(t,u)
#   rs[:,2] .> (rs[:,1] + rs[:,3])/2
# end



# sanity checks: does 2 st and 20 st yield a plausible distribution?
# 20 st does, but not 2, there are still periods of segregation
# something is wrong, probably about the input...

# the qualitivative feel of the first few figures is well matched
# I haven't checked the gamma distribution yet

# t,u = simulate(5,8,120)=

# plot([results[:,[a[1],a[2],a[3]]]])

# TODO: simulate and test bistable
