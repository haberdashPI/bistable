module Rankin2015

export BistableParams, Simulation, create_inputs, create_int_inputs,
  int_responses, responses, partial_collect, rlengths

using Parameters
using Base.Iterators: drop, take, partition
using Lazy

@with_kw struct BistableParams
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

struct Simulation
  t::Range{Float64}
  us::Base.Generator
  If::Matrix{Float64}
end

function create_inputs(DF,PR,sim_length;dt=1/1000,params=BistableParams())
  triplets = floor(Int,PR*sim_length/4)
  freqs = ones(Int,3*triplets)
  freqs[2:3:end] = 2
  times = 1/PR*ones(3*triplets)
  times[3:3:end] = 2/PR
  times[2:end] = cumsum(times[1:end-1])
  times[1] = 0.0

  create_input_helper(PR,0:dt:sim_length,dt,freqs,times,params)
end

function create_int_inputs(trial_length,DF,PR,sim_length;
                           dt=1/1000,params=BistableParams())

  triplets = floor(Int,PR*sim_length/4 - PR*sim_length/4/(trial_length+1))
  freqs = ones(Int,3*triplets)
  freqs[2:3:end] = 2

  times = 1/PR*ones(3*triplets)
  times[3:3:end] = 2/PR
  times[(3*trial_length):(3*trial_length):end] = 3/PR + 2/PR
  times[2:end] = cumsum(times[1:end-1])
  times[1] = 0.0

  create_input_helper(PR,0:dt:sim_length,dt,freqs,times,params)
end

function create_input_helper(PR,t,dt,freqs,times,params)
  const p = params

  function input(t::Float64)
    if t > 0.0
      (exp(2)/p.α_1^2) * t^2*exp(-2t / p.α_1) +
        (p.Λ*exp(2)/p.α_2^2) * t^2*exp(-2t / p.α_2)
    else 0.0
    end
  end

  If = zeros(length(t),2)
  max_t_len = floor(Int,10/dt*1/PR)
  for i in eachindex(freqs)
    start_t_i = floor(Int,times[i]/dt) + 1
    for t_i in start_t_i:(start_t_i + max_t_len)
      t_i <= length(t) || break
      If[t_i,freqs[i]] += input(t[t_i] - times[i])
    end
  end

  If
end

function Simulation(DF,PR,sim_length;dt=1/1000,params=BistableParams(),
                    inputs=create_inputs(DF,PR,sim_length,dt=dt,params=params))
  If = inputs

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
  xs = Matrix{Float64}(N,length(x))
  xs[1,:] = x
  i = 1
  for x in drop(itr,1)
    xs[i+=1,:] = x
    i < size(xs,1) || break
  end

  xs
end

const τ_response = 0.5

function response_fn(sim)
  let dt = sim.t[2] - sim.t[1],
    α = dt / (τ_response + dt)
    y = 0

    function response(u)
      y = α*u[r] + (1-α)y
      (y[1] + y[3])/2 > y[2]
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

struct Responses
  vals::Array{Bool}
  len::Float64
end

Base.mean(x::Responses) = mean(x.vals)

function responses(sim::Simulation)
  f = response_fn(sim)
  return Responses((f(u) for u in sim.us),sim.t[2] - sim.t[1])
end

function int_responses(trial_length,PR,sim::Simulation)
  dt = sim.t[2]- sim.t[1]
  unit = 1/PR * 4 * (trial_length+1)
  indices = @> (floor(Int,t/dt)+1 for t in 0:unit:maximum(sim.t)) drop(1)
  resp = response_fn(sim)

  state = start(sim.us)
  last = 1

  responses = map(indices) do i
    sum_x = sum(1:(i - last)) do j
      x,state = next(sim.us,state)
      resp(x)
    end
    last = i
    return sum_x > 0.5
  end

  Responses(responses,unit)
end

function rlengths(rs::Responses)
  oldval::Int = -1
  len::Int = 0
  lengths = Int[]
  for r in rs.vals
    if r*1 == oldval
      len += 1
    elseif len > 0
      push!(lengths,len)
      len = 1
      oldval = r*1
    else
      oldval = r*1
      len = 1
    end
  end
  push!(lengths,len)
  lengths .* rs.len
end

end
