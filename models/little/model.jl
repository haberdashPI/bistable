using HDF5
import Base: run

########################################
# util
function standardize!(x,dims)
  x .= x .- mean(x,dims)
  x .= x ./ std(x,dims)
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
struct Layer1
  w::Matrix{Float32}
  b::Vector{Float32}
  steps::Int
  inertia_n::Int
  inertia_threshold::Float32
end

function Layer1(filename)
  h5open(filename) do file
    Layer1(
      Float32.(read(file,"/layer1/W")),
      Float32.(squeeze(read(file,"/layer1/b")',2)),
      3,
      4,
      30
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

delta_t(m::Layer1) = m.steps * em.input_frame_len

sig(x,s,t) = 1 / (1 + exp(-4s*(x-t)))

function resemblance(m::Layer1,l1::Matrix{Float32},template::Matrix{Float32})
  template = reshapefor(template,m)
  standardize!(template,2)

  p = 0.25
  ytemp = (m.b' .+ template * m.w)[1,:]
  scale = norm(ytemp,p)
  mapslices(l1,2) do x_t
    # norm(x_t - ytemp,p) / scale
    sig(-norm(x_t - ytemp,p) / scale,10,-0.04)
    # 1 / (1 + exp(4(1 - norm(x_t .- ytemp) / ylen - 0.5)))
  end
end

function run(m::Layer1,x::Matrix{Float32})
  x = reshapefor(x,m)
  standardize!(x,2)

  y = @fastmath(m.b' .+ x * m.w)
  # # TODO: I could have multiple delta steps
  # # for a single frame to have more resolution in
  # # calculating the dynamics
  # y = similar(x,size(x,1),size(A,2))
  # y_old = zeros(x,size(m.w,2))
  # c = dt(m) / m.τ
  # for r in 1:floor(Int,size(x,1)/m.steps)
  #   y[r,:] = y_old
  #   y[r,:] .+= c .* (-y_old .+ m.b' .+ (x_i' * m.w))
  # end
  # y

  # TODO: do we keep the inertia? I think
  # I'm going to remove this when I introduce
  # some form of mutual inhibition
  inertia(m,y)
end

########################################
# layer 2

struct Layer2
  w::Vector{Matrix{Float32}}
  w_past::Vector{Matrix{Float32}}
  b::Vector{Vector{Float32}}
  steps::Int
end

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

ncomps(m::Layer2) = length(m.w)
steps(m::Layer2) = m.steps

function run(m::Layer2,x::Matrix{Float32})
  x = reshapefor(x,m)
  Σ = sum(1:ncomps(m)) do i
    y = x[2:end,:]*m.w[i]
    y .+= x[1:end-1,:]*m.w_past[i]
    y .+= m.b[i]'
    y
    @fast @views x[2:end,:]*m.w[i] .+ x[1:end-1,:]*m.w_past[i] .+ m.b[i]'
  end

  maxnorm!(Σ,2)
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
  spect_params::Vector{Float32}
  layer1::Layer1
  layer2::Vector{Layer2}
  layer3::Layer3
end

########################################
# model operations

const ntau = 4 #8
function Model(filename;spect_params=[10,8,-2,-1],thresh=0.9)
  layer1 = Layer1(filename)

  layer2 = Vector{Layer2}(ntau)
  for tau in 1:ntau
    info("Loading Layer2(tau=$tau)...")
    layer2[tau] = Layer2(filename,tau)
  end

  layer3 = Layer3(thresh)

  return Model(spect_params,layer1,layer2,layer3)
end

function run(model::Model,taus,x)
  x_l1 = run(model.layer1,Float32.(x))
  result = Vector{Matrix{Float32}}(maximum(taus))
  for tau in taus
    xl2_tau = run(model.layer2[tau],x_l1);
    result[tau] = run(model.layer3,xl2_tau)
  end
  result
end

frame_length(m::Model) = steps(m.layer1)

function respond(x,μ,σ²,N=1000)
  squeeze(mean(x .> μ + σ²*randn(size(x)...,N),ndims(x)+1),ndims(x)+1)
end
