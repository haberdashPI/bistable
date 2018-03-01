using ProgressMeter
using AxisArrays
using DSP
using HDF5
using Sounds

import DSP.Filters.freqs
import Sounds: Sound

export freqs, times, nfreqs, ntimes, delta_t, delta_f, Δt, Δf, frame_length,
  audiospect, Sound

const cochba = h5open(joinpath(@__DIR__,"..","data","cochba.h5")) do file
  read(file,"/real") + read(file,"/imag")*im
end

const fixed_fs = 8000

struct ASParams
  len::Int
  decay_tc::Float64
  nonlinear::Float64
  octave_shift::Float64
  fs::typeof(1.0Hz)
end

struct AuditorySpectrogram{T} <: ModelResult{T,2}
  val::AxisArray{T,2}
  params::ASParams
end
data(x::AuditorySpectrogram) = x.val
params(x::AuditorySpectrogram) = x.params
resultname(x::AuditorySpectrogram) = "Auditory Spectrogram"

nfreqs(as::ASParams) = size(cochba,2)-1
nfreqs(x) = length(freqs(x))
freqs(as::ASParams) = 440.0Hz * 2.0.^(((1:nfreqs(as)).-31)./24 .+ as.octave_shift)
freqs(as::ModelResult) = freqs(data(as))
freqs(as::AxisArray) = axisvalues(axes(as,Axis{:freq}))[1]

ntimes(x) = length(times(x))
times(as::ModelResult) = times(data(as))
times(as::AxisArray) = axisvalues(axes(as,Axis{:time}))[1]
times(p::ASParams,x::AbstractArray) = indices(x,1) .* Δt(p)

delta_t(x) = Δt(x)
delta_f(x) = Δf(x)

frame_length(params::ASParams) =
  round(Int,params.len * 2^(4+params.octave_shift))
Δt(params::ASParams) = uconvert(s,frame_length(params) / params.fs)
Δf(params::ASParams) = (1 / 24)*Hz
Sounds.samplerate(params::ASParams) = uconvert(Hz,params.fs)

frame_length(as::ModelResult) = frame_length(as.params)
Δt(as::ModelResult) = Δt(as.params)
Δf(as::ModelResult) = Δf(as.params)
Sounds.samplerate(x::ModelResult) = samplerate(params(x))

function ASParams(x;fs=samplerate(x),
                  len=10,decay_tc=8,nonlinear=-2,octave_shift=-1)
  @assert fs == fixed_fs*Hz "The only sample rate supported is $(fixed_fs)Hz"

  ASParams(len,decay_tc,nonlinear,octave_shift,fs)
end

function sigmoid(x::AbstractArray{T},fac::T) where T
  if fac > 0
    one(T) ./ (one(T) .+ exp.(.-x./fac))
  elseif fac == 0
    T.(x .> zero(T))
  elseif fac == -1
    max.(x,zero(T))
  elseif fac == -2
    T.(x)
    # TODO: implement halfregu
  else
    error("Non linear factor of $fac not supported")
  end
end


########################################
# auditory spectrogram
audiospect(x::AbstractArray;params...) = audiospect(x,ASParams(x;params...))

####################
# 'identity' conversions (something that's already basically a spectrogram)
function audiospect(x::AbstractArray,params::ASParams)
  if ndims(x) <= 2 && size(x,2) <= 2
    # the array probably represents a sound
    audiospect(Sound(x,rate=samplerate(params)))
  else
    # the array probably represents a spectrogram
    @assert(ndims(x) == 2,"Input to audiospect must be a sound or a "*
            "previously created spectrogram")
    @assert(size(x,2) == nfreqs(params),
            "Presumed input to be a spectrogram but the "*
            "number of columns do not match the number of frequency channels.")

    f = Axis{:freq}(freqs(params))
    t = Axis{:time}(times(params,x))
    AuditorySpectrogram(AxisArray(x,t,f),params)
  end
end

function audiospect(x::AxisArray{T,2} where T,params::ASParams)
  @assert(nfreqs(x) == nfreqs(params),
          "Frequency channels of array and parameters do not match")

  AuditorySpectrogram(x,params)
end

function audiospect(x::AuditorySpectrogram,params::ASParams)
  @assert(x.params == params,
          "Parameters of spectrogram and input parameters do not match")
  x
end

####################
# the actual computation of a spectrogram
audiospect(x::Sound,params::ASParams) =
  audiospect_helper(vec(Array(convert(Sound{fixed_fs,eltype(x),1},x))),params)

function audiospect_helper(x::Vector{T},params::ASParams,
                           internal_call=false) where {T}
  L, M = size(cochba)  # p_max = L - 2
  L_x = length(x)  # length of input
  frame_len  = frame_length(params)

  # decaying factor
  α = iszero(params.decay_tc) ? zero(T) :
    exp(-1/(params.decay_tc*2^(4+params.octave_shift)))

  # hair cell time constant in ms
  haircell_tc = 1//2
  β = exp(-1/(haircell_tc*2^(4+params.octave_shift)))

  N = ceil(Int,L_x / frame_len) # of frames
  if length(x) < N*frame_len
    append!(x,fill(zero(T),N*frame_len - length(x)))
  end
  v5 = fill(zero(T),N, M-1)
  y3_r = !internal_call ? fill(zero(T),0,0) : fill(zero(T),length(x),M-1)

  #######################################
  # last channel (highest frequency)
  #######################################

  p  = floor(Int,real(cochba[1, M]))
  B  = real(cochba[(0:p)+2, M])
  A  = imag(cochba[(0:p)+2, M])

  y1 = filt(PolynomialRatio(B,A),x)
  y2 = sigmoid(y1, params.nonlinear)

  # hair cell membrane (low-pass <= 4 kHz) ignored for LINEAR ionic channels
  if (params.nonlinear != -2) y2 = filt(PolynomialRatio([1.0],[1.0; -β]),y2) end
  y2_h = y2

  #######################################
  # All other channels
  #######################################
  for ch = (M-1):-1:1,

    #######################################
    # ANALYSIS: cochlear filterbank
    ########################################
    # (IIR) filter bank convolution ---> y1
    p  = floor(Int,real(cochba[1, ch]))  # order of ARMA filter
    B  = real(cochba[(0:p)+2, ch])  # moving average coefficients
    A  = imag(cochba[(0:p)+2, ch])  # autoregressive coefficients

    y1 = filt(B, A, x)
    ########################################
    # TRANSDUCTION: hair cells
    ########################################
    # Fluid cillia coupling (preemphasis) (ignored)

    # ionic channels (sigmoid function)
    y2 = sigmoid(y1, params.nonlinear)

    # hair cell membrane (low-pass <= 4 kHz) ---> y2 (ignored for linear)
    if (params.nonlinear != -2) y2 = filt(PolynomialRatio([1.0],[1.0; -β]),y2) end

    ########################################
    # REDUCTION: lateral inhibitory network
    ########################################
    # masked by higher (frequency) spatial response
    y3   = y2 - y2_h
    y2_h = y2
    if internal_call
      y3_r[:,ch] = y3
    end

    # spatial smoother ---> y3 (ignored)
    #y3s = y3 + y3_h
    #y3_h = y3

    # half-wave rectifier ---> y4
    y4 = max.(y3, 0)

    # temporal integration window ---> y5
    if !iszero(α)  # leaky integration
      y5 = filt(PolynomialRatio([1.0],[1.0; -α]),y4)
      v5[:, ch] = y5[frame_len*(1:N)]
    else    # short-term average
      if (frame_len == 1)
        v5[:, ch] = y4
      else
        v5[:, ch] = mean(reshape(y4, frame_len, N),1)'
      end
    end
  end

  if internal_call
    v5,y3_r
  else
    f = Axis{:freq}(freqs(params))
    t = Axis{:time}(times(params,v5))

    AuditorySpectrogram(AxisArray(v5,t,f),params)
  end
end

########################################
# inverse of auditory spectorgram
function Sounds.Sound(y_in::AuditorySpectrogram;max_iterations=typemax(Int),
                      target_error=0.05)
  @assert(max_iterations < typemax(Int) || target_error < Inf,
          "No stopping criterion specified (max_iterations or target_error).")
  params = y_in.params

  M = size(cochba,2)

  # expand y to include all frequencies
  y = zeros(eltype(y_in),size(y_in,1),M-1)
  f_ixs = minimum(freqs(y_in)) .<= freqs(y_in.params) .<= maximum(freqs(y_in))
  @assert sum(f_ixs) == size(y_in,2) "Unxpected frequency resolution."
  y[:,f_ixs] = y_in

  # generate initial guess
  x = inv_guess(params,y)

  # iteration setup
  ratios = similar(y)
  target_mean = mean(y)
  target_sum2 = sum(y.^2)

  min_err = Inf
  min_x = x

  if max_iterations < typemax(Int)
    prog = Progress(max_iterations,"Inverting Spectrogram: ")
  else
    prog = ProgressThresh(target_error,"Inverting Spectrogram: ")
  end

  for iteration in 1:max_iterations
    if min_err < target_error
      break
    end

    if params.nonlinear == 0
      x .-= mean(x)
      x ./= std(x)
    end

    y_hat,y3_hat = audiospect_helper(x,params,true)
    x = match_x(params,x,ratios,y,y_hat,y3_hat)

    y_hat .*= target_mean/mean(y_hat)
    err = sum((y_hat .- y).^2) ./ target_sum2

    if err < min_err
      min_x = x
      min_err = err
    elseif err-1 > min_err
      # restart
      x .= sign.(x) .+ rand(size(x))
      x .-= mean(x)
      x ./= std(x)
    end

    x .*= 1.01

    if max_iterations < typemax(Int)
      ProgressMeter.next!(prog;showvalues =
                          [(:error,string(100round(min_err,4),"%"))])
    else
      ProgressMeter.update!(prog,min_err)
    end
  end
  ProgressMeter.finish!(prog)

  Sound(min_x,rate=samplerate(params))
end

################################################################################
# private helper functions

function inv_guess(params::ASParams,y::AbstractMatrix)
  # ?? the initial guess only uses the first 48 channels
  f = ustrip.(uconvert.(Hz,freqs(params)))[1:48]
  steps = 1:frame_length(params)*size(y,1)
  indices = ceil.(Int,steps / frame_length(params))
  t = steps ./ ustrip(uconvert(Hz,params.fs))

  x = sum(cos.(2π.*f'.*t) .* view(y,indices,1:48),2)
  x .-= mean(x)
  x ./= std(x)

  squeeze(x,2)
end

function match_x(params::ASParams,x,ratios,y,y_hat,y3_hat)
  M = size(cochba,2)
  steps = 1:frame_length(params)*size(y,1)
  indices = ceil.(Int,steps / frame_length(params))

  map!(ratios,y_hat,y) do y_hat,y
    !iszero(y_hat) ? y ./ y_hat : !iszero(y) ? 2 : 1
  end

  x .= 0
  for ch in 1:M-1
    p = floor(Int,real(cochba[1, ch]))  # order of ARMA filter
    ch_norm  = imag(cochba[1, M])
    B = real(cochba[(0:p)+2, ch])  # moving average coefficients
    A = imag(cochba[(0:p)+2, ch])  # autoregressive coefficients

    if params.nonlinear == -2
      y1 = y3_hat[:,ch].*view(ratios,indices,ch)
    else
      y1 = y3_hat[:,ch]
      posi = find(y1 .>= 0)
      y1[posi] .*= view(ratios,indices[posi],ch)

      negi = find(y1 .< 0)
      y1[negi] .*= maximum(y1[posi]) / -minimum(y1[negi])
    end

    x .+= reverse(filt(B,A,reverse(y1))) / ch_norm
  end

  x
end

function freq_ticks(as)
  a = minimum(freqs(as))
  b = maximum(freqs(as))
  step = 0.25

  helper(step) = round.(filter(f -> a <= f*Hz <= b,1000*2.0.^(-3:step:2)),-1)
  fbreaks = helper(step)
  while length(fbreaks) > 7
    fbreaks = helper(step *= 2)
  end

  fs = ustrip(uconvert.(Hz,freqs(as)))

  findices = mapslices(abs.(fbreaks .- fs'),2) do row
    _, i = findmin(row)
    i
  end

  fbreaks,findices
end
