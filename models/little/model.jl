using HDF5
using Unitful: s
using Parameters
import Base: run

# TODO: what if we include the sigmoid operations (the reall model) at each
# layer instead of the weird, maxnorm thing that is happening right now?

include("audio_spect.jl")

const Seconds = typeof(1f0*s)
########################################
# util
function standardize!(x,dims...)
  x .= x .- mean(x,dims...)
  x .= x ./ max.(1e-16,std(x,dims...))
end

function maxnorm!(x,dims)
  x ./= max.(1e-16,abs.(maximum(x,dims)))
end

function reshapefor(x,model)
  n = floor(Int,size(x,1)/steps(model))
  reshape(x[1:(steps(model)*n),:]',:,n)'
end

@with_kw struct LayerParams
  τ_a::Seconds = 1.5s
  c_a::Float32 = 1
  c_mi::Float32 = 0.6
  use_sig::Bool = false
end

sig(x) = 1 / (1 + exp(-x))

function adapt(fn,layer,x,indices=1:size(x,1))
  p = layer.params
  y_t = fill(zero(x[1]),size(layer.w,2))
  c_a = delta_t(layer) / p.τ_a
  c_mi = p.c_mi / length(y_t)
  a = zeros(y_t)

  y = similar(x,size(x,1),size(layer.w,2))

  for t in indices
    y_t .= fn(x,t) .-
      p.c_a.*a .- # adaptation
      (y_t .- sum(y_t)).*c_mi # mutual inhibitoin
    y[t,:] = p.use_sig ? y_t .= sig.(y_t) : y_t

    a .+= c_a.*(y_t .- a)
  end

  y
end

########################################
# layer 1
struct Layer1
  w::Matrix{Float32}
  b::Vector{Float32}
  bv::Vector{Float32}

  steps::Int
  params::LayerParams
  input_frame_len::Seconds

  inertia_n::Int
  inertia_threshold::Float32
end

Base.show(io::IO,x::Layer1) =
  write(io,"Layer1(τ_a=$(x.params.τ_a),"*
        "c_a=$(x.params.c_a),inertia_n=$(x.inertia_n),"*
        "intertia_θ=$(x.inertia_threshold))")

function Layer1(filename,audio_spect,params=LayerParams())
  h5open(filename) do file
    Layer1(
      Float32.(read(file,"/layer1/W")),
      Float32.(squeeze(read(file,"/layer1/b")',2)),
      Float32.(squeeze(read(file,"/layer1/bv")',2)),
      3,params,
      1/audio_spect.fs * frame_length(audio_spect) * s,4,30
    )
  end
end
steps(m::Layer1) = m.steps;

function inertia(m::Layer1,out::Matrix{Float32})
  start = view(out,1:m.inertia_n,:)
  map!(start,start) do x
    x > m.inertia_threshold ? 5.0x : x/5.0
  end

  out
end

unit_ordering_by(m::Layer1,from,to) =
  unit_ordering_by(m,Float32.(from),Float32.(to))
function unit_ordering_by(m::Layer1,
                          from::Matrix{Float32},
                          to::Matrix{Float32})
  function resp(x)
    x = reshapefor(x,m)
    standardize!(x,2)
    m.b' .+ x * m.w
  end

  levels = abs.(resp(from) .- resp(to))

  sort(indices(m.w,2),by = i -> -levels[1,i])
end

run(m::Layer1,x::Matrix;params...)= run(m,Float32.(x);params...)

delta_t(m::Layer1) = m.steps * m.input_frame_len

function run(m::Layer1,x::Matrix{Float32};use_inertia=true)
  x = reshapefor(x,m)
  standardize!(x,2)
  p = m.params

  y = adapt(m,x) do x,t
    (x[t,:]' * m.w)' .+ m.b
  end

  use_inertia ? inertia(m,y) : y
end

# get a visible response for each hidden unit
function display(m::Layer1)
  sig.(m.w .+ m.bv)
end

########################################
# layer 2

struct Layer2
  w::Array{Float32,3}
  w_past::Array{Float32,3}
  b::Array{Float32,2}
  bv::Array{Float32,2}
  params::LayerParams
  steps::Int
  layer1_len::Seconds
  function Layer2(w,w_past,b,bv,params,steps,layer1_len)
    @assert(size(w,3) == size(w_past,3) == size(b,2) == size(bv,2),
            "w, w_past and b lengths must be equal")
    new(w,w_past,b,bv,params,steps,layer1_len)
  end
end
Base.show(io::IO,x::Layer2) =
  !get(io,:compact,false) ? write(io,"Layer2(τ=$(x.steps))") :
  write(io,"L2(τ=$(x.steps))")

function Layer2(filename::String,layer1::Layer1,steps::Int,params::LayerParams)
  h5open(filename) do file
    Layer2(
      Float32.(read(file,"/layer2/tau$steps/W")),
      Float32.(read(file,"/layer2/tau$steps/Wpast")),
      Float32.(read(file,"/layer2/tau$steps/b")),
      Float32.(read(file,"/layer2/tau$steps/bv")),
      params,
      steps,
      delta_t(layer1)
    )
  end
end

steps(m::Layer2) = m.steps
delta_t(m::Layer2) = m.layer1_len * m.steps

function run(m::Layer2,x::Matrix{Float32})
  x = reshapefor(x,m)
  p = m.params
  y_sum = sum(1:size(m.w,3)) do i
    y = adapt(m,x,2:size(x,1)) do x,t
      (view(x,t,:)'*m.w[:,:,i])' .+
        (view(x,t-1,:)'*m.w_past[:,:,i])' .+
        m.b[:,i]
    end
    @views y[2:end,:]
  end
  # NOTE: as written this just assumes z_k values are uniform.
  # this is a bug, but I want to change things slowly from the original
  # deb model. (It's also not clear how adaptation would apply to
  # during the calculation of z).
  p.use_sig ? y_sum / size(m.w,3) : maxnorm!(y_sum,2)
end

# get a visible response for each hidden unit at specified set of
# time steps for a given mixture
function display(m::Layer2,x::Matrix{Float32},i)
  p = resize((x[1:end-1]'*m.w_past[:,:,i]),size(m.w_past,1),1,size(x,1)-1)
  m.w[:,:,i] .+ p .+ m.bv
end

########################################
# layer 3

struct Layer3
  thresh::Float32
  params::LayerParams
  delta_t::Seconds
end

# function Layer3(thresh::Float32,params::LayerParams,layer2::Layer2)
#   Layer3(thresh,params,delta_t(layer2))
# end

delta_t(m::Layer3) = m.delta_t

function run(m::Layer3,x::Matrix{Float32})
  p = m.params
  resp = zeros(Float32,size(x))
  C = zeros(Float32,size(x,2),size(x,2))
  A = zeros(C)
  c_mi = p.c_mi / length(C)
  c_a = delta_t(m) / p.τ_a
  for i in 1:size(x,1)
    k = x[i,:] .* sign.(x[i,:].-m.thresh)
    # TODO: create the adaptation and MI dynamics here
    C .+= k .* k' .-
      p.c_a .* A .-
      (C .- sum(C)).*c_mi
      # for j in 1:size(C,1) C[j,j] = 0.0 end

    A .+= c_a.*(C .- A)

    resp[i,:] = 2maxnorm!(C*x[i,:],1)
  end

  resp
end

########################################
# model defintion

struct Model
  spect::AuditorySpectrogram
  layer1::Layer1
  layer2::Vector{Layer2}
  layer3::Vector{Layer3}
end

########################################
# model operations

const ntau = 4 #8
function Model(dir;
               l1params = LayerParams(),
               l2params = LayerParams(),
               l3params = LayerParams(),
               layer3_thresh=0.9,
               spect_len=10,spect_decay_tc=8,
               spect_nonlinear=-2,spect_octave_shift=-1)

  audio_spect = AuditorySpectrogram(
    joinpath(dir,"cochba.h5"),
    len=spect_len,
    decay_tc=spect_decay_tc,
    nonlinear=spect_nonlinear,
    octave_shift=spect_octave_shift
  )

  filename = joinpath(dir,"model.h5")
  layer1 = Layer1(filename,audio_spect,l1params)

  layer2 = Vector{Layer2}(ntau)
  layer3 = Vector{Layer3}(ntau)
  for tau in 1:ntau
    info("Loading Layer2(tau=$tau)...")
    layer2[tau] = Layer2(filename,layer1,tau,l2params)
    layer3[tau] = Layer3(layer3_thresh,l3params,delta_t(layer2[tau]))
  end

  return Model(audio_spect,layer1,layer2,layer3)
end

# avoid redefining this constant when the script is re-evaluated
@static if !isdefined(:spect_cache)
  const spect_cache = Dict{Vector,Matrix}()
end

function run_spect(model::Model,x::Vector)
  get!(spect_cache,x) do
    run(model.spect,x)
  end
end

run(model::Model,taus,x::AbstractArray;keys...) =
  run(model,taus,run_spect(model,Float64.(x));keys...)
function run(model::Model,taus,x::Matrix;upto=3)
  if upto == 0
    return x
  end

  l1 = run(model.layer1,Float32.(x))
  if upto == 1
    return l1
  end
  l2 = Vector{Matrix{Float32}}(maximum(taus))
  l3 = Vector{Matrix{Float32}}(maximum(taus))
  for tau in taus
    l2[tau] = run(model.layer2[tau],l1)
    if upto > 2
      l3[tau] = run(model.layer3[tau],l2[tau])
    end
  end
  if upto == 2
    l2
  else
    l3
  end
end

layer1_ordering_by(m::Model,x::Vector,y::Vector) =
  unit_ordering_by(m.layer1,Float32.(run_spect(m,x)),Float32.(run_spect(m,y)))

frame_length(m::Model) = s/m.spect.fs * steps(m.layer1) * frame_length(m.spect)

function respond(x,μ,σ²,N=1000)
  squeeze(mean(x .> μ + σ²*randn(size(x)...,N),ndims(x)+1),ndims(x)+1)
end
