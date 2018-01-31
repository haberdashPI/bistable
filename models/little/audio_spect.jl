using ProgressMeter
using RCall
using DataFrames
using DSP
using HDF5
import Base: run
using Sounds
using RecipesBase

import DSP.Filters.freqs

include("plots.jl")

struct AuditorySpectrogram
  cochba::Matrix{Complex128}
  len::Int
  decay_tc::Float64
  nonlinear::Float64
  octave_shift::Float64
  fs::Int
  min_freq::Hertz{Float64}
  max_freq::Hertz{Float64}
end

SampledSignals.nchannels(as::AuditorySpectrogram) = size(as.cochba,2)-1
channels_computed(s::AuditorySpectrogram) =
  find(f -> s.min_freq <= f <= s.max_freq,all_freqs(s))

Base.show(io::IO,x::AuditorySpectrogram) =
  write(io,"AuditorySpectrogram(len=$(x.len),decay_tc=$(x.decay_tc),"*
        "nonlinear=$(x.nonlinear),octave_shift=$(x.octave_shift))")

function AuditorySpectrogram(filename::String;
                             fs=ustrip(samplerate()),
                             len=10,decay_tc=8,nonlinear=-2,octave_shift=-1,
                             min_freq = -Inf*Hz,max_freq = Inf*Hz)
  min_freq = convert(Hertz{Float64},min_freq)
  max_freq = convert(Hertz{Float64},max_freq)

  @assert fs == 8000 "The only sample rate supported is 8000 Hz"
  h5open(filename) do file
    r = read(file,"/real")
    i = read(file,"/imag")
    AuditorySpectrogram(r + i*im,len,decay_tc,nonlinear,octave_shift,fs,
                        min_freq,max_freq)
  end
end

all_freqs(as::AuditorySpectrogram) =
  440.0Hz * 2.0.^(((1:nchannels(as)).-31)./24 .+ as.octave_shift)

freqs(as::AuditorySpectrogram) = freqs(as,1:length(channels_computed(as)))
freqs(as::AuditorySpectrogram,data::AbstractMatrix) = freqs(as,indices(data,2))
freqs(as::AuditorySpectrogram,data::AbstractArray{T,4}) where T =
    freqs(as,indices(data,2))
freqs(as::AuditorySpectrogram,ixs) =
  filter(f -> as.min_freq <= f <= as.max_freq,all_freqs(as))[ixs]

times(as::AuditorySpectrogram,data::AbstractMatrix) =
  indices(data,1) ./ as.fs .* frame_length(as) * s

@recipe function plot(as::AuditorySpectrogram,data::Matrix)
  @series begin
    seriestype := :heatmap
    fillcolor --> :reds
    xlabel --> "Time (s)"
    ylabel --> "Frequency (Hz)"
    (ustrip(times(as,data)),ustrip(freqs(as,data)),data')
  end
end

function freq_ticks(as::AuditorySpectrogram,x)
  a = ustrip(as.min_freq)
  b = ustrip(as.max_freq)
  step = 0.25

  helper(x,step) = round.(filter(f -> a <= f <= b,1000*2.0.^(-3:step:2)),-1)
  fbreaks = helper(x,step)
  while length(fbreaks) > 7
    fbreaks = helper(x,(step *= 2))
  end

  fs = ustrip(freqs(as,x))

  findices = mapslices(abs.(fbreaks .- fs'),2) do row
    _, i = findmin(row)
    i
  end

  fbreaks,findices
end

rplot(as::AuditorySpectrogram,data::Sound) = rplot(as,as(data))
rplot(as::AuditorySpectrogram,data::AbstractVector) = rplot(as,as(data))
function rplot(as::AuditorySpectrogram,data::Matrix)
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(as,data)[at(ixs,1)])),
                 freq_bin = vec(at(ixs,2)))
  fbreaks,findices = freq_ticks(as,data)
  p = raster_plot(df,value=:response,x=:time,y=:freq_bin)
R"""

  library(ggplot2)

  $p +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (Hz)') + xlab('Time (s)')

"""
end

function sigmoid(x,fac)
  if fac > 0
    1 ./ (1 .+ exp.(-x./fac))
  elseif fac == 0
    x .> 0.0
  elseif fac == -1
    max.(x,0.0)
  elseif fac == -2
    x
    # TODO: implement halfregu
  else
    error("Non linear factor of $fac not supported")
  end
end

frame_length(s::AuditorySpectrogram) = round(Int,s.len * 2^(4+s.octave_shift))
Δt(as::AuditorySpectrogram) = s * frame_length(as) / as.fs
Δf(as::AuditorySpectrogram) = 1 / 24

(s::AuditorySpectrogram)(x::Sound) = s(Array(convert(Sound{s.fs,Float64,1},x)))

function (s::AuditorySpectrogram)(x::Vector{T},internal_call=false) where T
  L, M = size(s.cochba)	# p_max = L - 2
  L_x = length(x)	# length of input
  frame_len	= frame_length(s)

  # decaying factor
  α = iszero(s.decay_tc) ? zero(T) : exp(-1/(s.decay_tc*2^(4+s.octave_shift)))

  # hair cell time constant in ms
  haircell_tc = 1//2
  β = exp(-1/(haircell_tc*2^(4+s.octave_shift)))

  N = ceil(Int,L_x / frame_len) # of frames
  if length(x) < N*frame_len
    append!(x,fill(zero(T),N*frame_len - length(x)))
  end
  v5 = fill(zero(T),N, M-1)
  y3_r = !internal_call ? fill(zero(T),0,0) : fill(zero(T),length(x),M-1)

  #######################################
  # last channel (highest frequency)
  #######################################

	p	= floor(Int,real(s.cochba[1, M]))
	B	= real(s.cochba[(0:p)+2, M])
	A	= imag(s.cochba[(0:p)+2, M])

  y1	= filt(PolynomialRatio(B,A),x)
  y2	= sigmoid(y1, s.nonlinear)

  # hair cell membrane (low-pass <= 4 kHz) ignored for LINEAR ionic channels
  if (s.nonlinear != -2) y2 = filt(PolynomialRatio([1.0],[1.0; -β]),y2) end
  y2_h = y2

  #######################################
  # All other channels
  #######################################
  for ch = (M-1):-1:1,

	  #######################################
	  # ANALYSIS: cochlear filterbank
	  ########################################
	  # (IIR) filter bank convolution ---> y1
	  p  = floor(Int,real(s.cochba[1, ch]))	# order of ARMA filter
	  B  = real(s.cochba[(0:p)+2, ch])	# moving average coefficients
	  A  = imag(s.cochba[(0:p)+2, ch])	# autoregressive coefficients

	  y1 = filt(B, A, x)
	  ########################################
	  # TRANSDUCTION: hair cells
	  ########################################
	  # Fluid cillia coupling (preemphasis) (ignored)

	  # ionic channels (sigmoid function)
	  y2 = sigmoid(y1, s.nonlinear)

	  # hair cell membrane (low-pass <= 4 kHz) ---> y2 (ignored for linear)
	  if (s.nonlinear != -2) y2 = filt(PolynomialRatio([1.0],[1.0; -β]),y2) end

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
	  if !iszero(α)	# leaky integration
		  y5 = filt(PolynomialRatio([1.0],[1.0; -α]),y4)
		  v5[:, ch] = y5[frame_len*(1:N)]
	  else		# short-term average
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
    v5[:,channels_computed(s)]
  end
end

function inv_guess(spect::AuditorySpectrogram,y::AbstractMatrix)
  # ?? the initial guess only uses the first 48 channels
  f = ustrip(all_freqs(spect))[1:48]
  steps = 1:frame_length(spect)*size(y,1)
  indices = ceil.(Int,steps / frame_length(spect))
  t = steps ./ spect.fs

  x = sum(cos.(2π.*f'.*t) .* view(y,indices,1:48),2)
  x .-= mean(x)
  x ./= std(x)

  squeeze(x,2)
end

function match_x(spect::AuditorySpectrogram,x,ratios,y,y_hat,y3_hat)
  M = size(spect.cochba,2)
  steps = 1:frame_length(spect)*size(y,1)
  indices = ceil.(Int,steps / frame_length(spect))

  map!(ratios,y_hat,y) do y_hat,y
    !iszero(y_hat) ? y ./ y_hat : !iszero(y) ? 2 : 1
  end

  x .= 0
  for ch in 1:M-1
	  p = floor(Int,real(spect.cochba[1, ch]))	# order of ARMA filter
		ch_norm	= imag(spect.cochba[1, M])
	  B = real(spect.cochba[(0:p)+2, ch])	# moving average coefficients
	  A = imag(spect.cochba[(0:p)+2, ch])	# autoregressive coefficients

    if spect.nonlinear == -2
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

function Base.inv(spect::AuditorySpectrogram,y_in::AbstractMatrix;
                  max_iterations=typemax(Int),max_error=0.05,usematlab=false)
  @assert(max_iterations < typemax(Int) || max_error < Inf,
          "No stopping criterion specified (max_iterations or max_error).")

  y = y_in
  M = size(spect.cochba,2)

  # expand y to include all frequencies
  y = similar(y_in,size(y,1),M-1)
  y[:,channels_computed(spect)] = y_in

  # generate initial guess
  x = inv_guess(spect,y)

  # iteration setup
  ratios = similar(y)
  target_mean = mean(y)
  target_sum2 = sum(y.^2)

  min_err = Inf
  min_x = x

  if max_iterations < typemax(Int)
    prog = Progress(max_iterations,"Inverting Spectrogram: ")
  else
    prog = ProgressThresh(max_error,"Inverting Spectrogram: ")
  end

  for iteration in 1:max_iterations
    if min_err < max_error
      break
    end

    if spect.nonlinear == 0
      x .-= mean(x)
      x ./= std(x)
    end

    y_hat,y3_hat = spect(x,true)
    x = match_x(spect,x,ratios,y,y_hat,y3_hat)

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

  min_x
end
