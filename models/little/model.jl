using HDF5
using Unitful: s
import Base: run

include("audio_spect.jl")
########################################
# util
function standardize!(x,dims...)
  x .= x .- mean(x,dims...)
  x .= x ./ std(x,dims...)
end

function maxnorm!(x,dims)
  x ./= abs.(maximum(x,dims))
end

function reshapefor(x,model)
  n = floor(Int,size(x,1)/steps(model))
  reshape(x[1:(steps(model)*n),:]',:,n)'
end

########################################
# layer 1

const Seconds = typeof(1f0*s)

struct Layer1
  w::Matrix{Float32}
  b::Vector{Float32}

  steps::Int
  N_a::Int
  τ_a::Seconds
  c_a::Float32

  input_frame_len::Seconds

  inertia_n::Int
  inertia_threshold::Float32
end
Base.show(io::IO,x::Layer1) =
  write(io,"Layer1(N_a=$(x.N_a),τ_a=$(x.τ_a),c_a=$(x.c_a),"*
        "inertia_n=$(x.inertia_n),intertia_θ=$(x.inertia_threshold))")

function Layer1(filename,audio_spect;
                N_a=5,τ_a=1.5s,c_a=0.3)
  h5open(filename) do file
    Layer1(
      Float32.(read(file,"/layer1/W")),
      Float32.(squeeze(read(file,"/layer1/b")',2)),
      3,N_a,τ_a,c_a,
      1/audio_spect.fs * frame_length(audio_spect) * s,4,30
    )
  end
end
steps(m::Layer1) = m.steps;

function inertia(m::Layer1,out::Matrix{Float32})
  start = view(out,1:m.inertia_n,:)
  helper(x) = x > m.inertia_threshold ? 5.0x : x/5.0
  start .= helper.(start)
  return out
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

function run(m::Layer1,x::Matrix{Float32};
             use_inertia=true,return_adaptation=false)
  N_t = floor(Int,size(x,1)/steps(m))
  y = similar(x,N_t,size(m.w,2))
  adapt = zeros(Float32,N_t*m.N_a+1,size(m.w,2))
  adapt_index = 1
  c = delta_t(m) / m.N_a / m.τ_a

  old = similar(x,size(m.w,2))
  for t in 1:N_t
    indices = m.steps*(t-1)+1 : m.steps*t

    # TODO: if we redefine the weight ordering, we could reorganize
    # weights to make this transpose unncessary
    x_i = standardize!(x[indices,:]'[:])
    y[t,:] .= m.b .+ (x_i' * m.w)' .- m.c_a.*adapt[adapt_index,:]

    old .= view(adapt,adapt_index,:)
    for i in 1:m.N_a
      adapt_index += 1
      old .= adapt[adapt_index,:] .= old .+ c.*(y[t,:] .- old)
    end
  end

  y = use_inertia ? inertia(m,y) : y

  return_adaptation ? (y,adapt) : y
end

########################################
# layer 2

struct Layer2
  w::Vector{Matrix{Float32}}
  w_past::Vector{Matrix{Float32}}
  b::Vector{Vector{Float32}}
  steps::Int
  function Layer2(w,w_past,b,steps)
    @assert(length(w) == length(w_past) == length(b),
            "w, w_past and b lengths must be equal")
    new(w,w_past,b,steps)
  end
end
Base.show(io::IO,x::Layer2) =
  !get(io,:compact,false) ? write(io,"Layer2(τ=$(x.steps))") :
  write(io,"L2(τ=$(x.steps))")

function Layer2(filename::String,i::Int)
  function helper(x::Array{T,3}) where T
    collect(Float32.(x[:,:,i]) for i in 1:size(x,3))
  end

  function helper(x::Array{T,2}) where T
    collect(Float32.(x[:,i]) for i in 1:size(x,2))
  end

  h5open(filename) do file
    Layer2(
      helper(read(file,"/layer2/tau$i/W")),
      helper(read(file,"/layer2/tau$i/Wpast")),
      helper(read(file,"/layer2/tau$i/b")),
      i
    )
  end
end

steps(m::Layer2) = m.steps

function run(m::Layer2,x::Matrix{Float32})
  x = reshapefor(x,m)
  y = sum(1:length(m.w)) do i
    @fastmath @views x[2:end,:]*m.w[i] .+ x[1:end-1,:]*m.w_past[i] .+ m.b[i]'
  end
  maxnorm!(y,2)
end

########################################
# layer 3

struct Layer3
  thresh::Float32
end

function run(m::Layer3,x::Matrix{Float32})
  resp = zeros(Float32,size(x))
  C = zeros(Float32,size(x,2),size(x,2))
  @fastmath @views for i in 1:size(x,1)
    k = x[i,:] .* sign.(x[i,:].-m.thresh)
    C .+= k .* k'
    # for j in 1:size(C,1) C[j,j] = 0.0 end

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
  layer3::Layer3
end

########################################
# model operations

const ntau = 4 #8
function Model(dir;
               l1_N_a=1,
               l1_τ_a=1.5s,
               l1_c_a=1f0,
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
  layer1 = Layer1(filename,audio_spect,
                  N_a=l1_N_a,
                  τ_a=l1_τ_a,
                  c_a=l1_c_a)

  layer2 = Vector{Layer2}(ntau)
  for tau in 1:ntau
    info("Loading Layer2(tau=$tau)...")
    layer2[tau] = Layer2(filename,tau)
  end

  layer3 = Layer3(layer3_thresh)

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

run(model::Model,taus,x::Vector) = run(model,taus,run_spect(model,x))
function run(model::Model,taus,x::Matrix)
  x_l1 = run(model.layer1,Float32.(x))
  result = Vector{Matrix{Float32}}(maximum(taus))
  for tau in taus
    xl2_tau = run(model.layer2[tau],x_l1);
    result[tau] = run(model.layer3,xl2_tau)
  end
  result
end

layer1_ordering_by(m::Model,x::Vector,y::Vector) =
  unit_ordering_by(m.layer1,Float32.(run_spect(m,x)),Float32.(run_spect(m,y)))

frame_length(m::Model) = s/m.spect.fs * steps(m.layer1) * frame_length(m.spect)

function respond(x,μ,σ²,N=1000)
  squeeze(mean(x .> μ + σ²*randn(size(x)...,N),ndims(x)+1),ndims(x)+1)
end
