using RCall
using DataFrames
using MATLAB
using DSP
using HDF5
import Base: run
using TimedSound
using RecipesBase

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

Base.show(io::IO,x::AuditorySpectrogram) =
  write(io,"AuditorySpectrogram(len=$(x.len),decay_tc=$(x.decay_tc),"*
        "nonlinear=$(x.nonlinear),octave_shift=$(x.octave_shift))")

base_min_freq(shift) = Hz*1000.*2.^(1 / 24 - 2.5)
base_max_freq(shift) = Hz*1000.*2.^(128 / 24 - 2.5)

function AuditorySpectrogram(filename::String;
                             fs=ustrip(TimedSound.samplerate()),
                             len=10,decay_tc=8,nonlinear=-2,octave_shift=-1,
                             min_freq = base_min_freq(octave_shift),
                             max_freq = base_max_freq(octave_shift))
  mat"loadload;"
  min_freq = convert(Hertz{Float64},min_freq)
  max_freq = convert(Hertz{Float64},max_freq)

  @assert(min_freq >= base_min_freq(octave_shift),
          "$min_freq is lower than the minimum allowed ($(base_min_freq(octave_shift)))")
  @assert(max_freq <= base_max_freq(octave_shift),
          "$max_freq is higher than the maximum allowed ($(base_max_freq(octave_shift)))")

  @assert fs == 8000 "The only sample rate supported is 8000 Hz"
  h5open(filename) do file
    r = read(file,"/real")
    i = read(file,"/imag")
    AuditorySpectrogram(r + i*im,len,decay_tc,nonlinear,octave_shift,fs,
                        min_freq,max_freq)
  end
end

all_freqs(as::AuditorySpectrogram) =
  440.0Hz * 2.0.^(((1:size(as.cochba,2)).-31)./24 .+ as.octave_shift)

freqs(as::AuditorySpectrogram) = freqs(as,1:size(as.cochba,2)-1)
freqs(as::AuditorySpectrogram,data::AbstractMatrix) = freqs(as,indices(data,2))
freqs(as::AuditorySpectrogram,ixs) =
  filter(f -> as.min_freq <= f <= as.max_freq,all_freqs(as))[ixs]

times(as::AuditorySpectrogram,data::AbstractMatrix) =
  indices(data,1) ./ as.fs .* frame_length(as) * s

@recipe function plot(as::AuditorySpectrogram,data::Matrix)
  @series begin
    seriestype := :heatmap
    (times(as,data),freqs(as,data),data')
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

rplot(as::AuditorySpectrogram,data::TimedSound.Sound) = rplot(as,as(data))
rplot(as::AuditorySpectrogram,data::AbstractVector) = rplot(as,as(data))
function rplot(as::AuditorySpectrogram,data::Matrix)
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(as,data)[at(ixs,1)])),
                 freq_bin = vec(at(ixs,2)))
  fbreaks,findices = freq_ticks(as,data)
R"""

  library(ggplot2)

  ggplot($df,aes(x=time,y=freq_bin,fill=response)) +
    geom_raster() +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (Hz)') + xlab('Time (s)') +
    scale_fill_distiller(palette='Reds',direction=1)

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

(s::AuditorySpectrogram)(x::TimedSound.Sound{8000,T}) where T =
  s(Float64.(x[:,:left]))
(s::AuditorySpectrogram)(x::TimedSound.Sound{R}) where R =
  error("sound must be 8kHz mono.")

function (s::AuditorySpectrogram)(x::Vector{T}) where T
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

  ixs = find(f -> s.min_freq <= f <= s.max_freq,all_freqs(s))
  v5[:,ixs]
end

function Base.inv(spect::AuditorySpectrogram,y::AbstractMatrix;iterations=10,
                  usematlab=true)
  paras = [spect.len, spect.decay_tc, spect.nonlinear, spect.octave_shift,
           iterations, 0, 0]
  my = mxarray(Float64.(y))
  guess = mat"aud2wavi($my,$paras)"
  mat"aud2wav($my,$guess,$paras)"

  # for julia implementation
  # 1. run spect forward with current guess

  # 2. for each channel, create new scaling values by comparing to target
  #    spectrum, and apply them

  # 3. apply any necessary hair cell modification

  # 4. reverse filter and apply result to y_cum (the eventually output)

  # 5. rescale output and calculate error ('marked previous performance')
end
