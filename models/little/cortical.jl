using RecipesBase
using DataFrames
using RCall
include("audio_spect.jl")

abstract type CorticalRepresentation end

struct CorticalModel <: CorticalRepresentation
  aspect::AuditorySpectrogram
  rates::Vector{Float64}
  scales::Vector{Float64}
  bandonly::Bool
end

scales(cm::CorticalModel) = cm.scales
rates(cm::CorticalModel) = [.-cm.rates; cm.rates]
freqs(cm::CorticalModel,data::Array{T,4}) where T =
    @views freqs(cm.aspect,data[:,1,1,:])
times(cm::CorticalModel,data::Array{T,4}) where T =
    @views times(cm.aspect,data[:,1,1,:])

function CorticalModel(aspect::AuditorySpectrogram;
                       rates=2.^(1:0.5:5),scales=2.^(-2:0.5:3),
                       bandonly=true)
  CorticalModel(aspect,rates,scales,true)
end

Δt(c::CorticalModel) = Δt(c.aspect)

@recipe function plot(cm::CorticalModel,data::Array{T,4}) where T
  N_r, N_s = length(cm.rates), length(cm.scales)
  layout := (N_r, N_s)
  for (i,x) in enumerate(CartesianRange((N_r,N_s)))
    r,s = x.I

    x = @views times(cm.aspect,data[:,1,1,:])
    y = @views freqs(cm.aspect,data[:,1,1,:])

    @series begin
      seriestype := :heatmap
      subplot := i
      (x,y,data[:,r,s,:]')
    end
  end
end

const spect_rate = 24

pad(x,lens) = pad(x,lens...)
function pad(x,lens::T...) where T <: Number
  @assert all(size(x) .<= lens)
  y = zeros(eltype(x),lens)
  y[indices(x)...] = x
  y
end

# transforms a bandpass frequency response into either a high or low pass
# response (or leaves it untouched)
function askind(H,len,maxi,kind)
  if kind == :low
    old_sum = sum(H)
    H[1:maxi-1] = 1
    H .= H ./ sum(H) .* old_sum
    H
  elseif kind == :high
    old_sum = sum(H)
    H[maxi+1:len] = 1
    H .= H ./ sum(H) .* old_sum
    H
  else
    H
  end
end

# create the frequency-scale filter (filter along spectral axis)
function scale_filter(scale,len,ts,kind)
  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)
  _,maxi = findmax(H)

  askind(H,len,maxi,kind)
end

# create the temporal-rate filter (filter along temporal axis)
function rate_filter(rate,len,spect_len,kind)
  t = (0:len-1)./spect_len .* rate
  h = @. rate * t^3 * exp(-4t) * cos(2π*t)
  h .-= mean(h)

  H0 = view(fft(pad(h,2len)),1:len)
  A = angle.(H0)
  H = abs.(H0)

  _, maxi = findmax(H)
  H ./= H[maxi]
  H = askind(H,len,maxi,kind)

  H .* exp.(A*im)
end

function (cm::CorticalModel)(s::Matrix;rates=cm.rates,scales=cm.scales)
  N_t, N_f = size(s)
  N_r, N_s = length(rates), length(scales)

  # spatial-frequency represention of the spectrogram (i.e. 2D fft of s).

  # implementation note: we transpose s because rfft halves the first dimension
  # and we want frequency not time to be the halved dimension (since the rate
  # filters, over time, need both negative and positive direction and that
  # invovles manipulating negative frequencies)
  S = rfft(pad(s',2.*reverse(size(s))))
  st_ifft = plan_irfft(S,2size(s,2))

  cr = zeros(eltype(s), N_t, 2N_r, N_s, N_f)
  for (ti,rate) in enumerate(rates)
	  # rate filtering
	  HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                         cm.bandonly ? :band :
                         ti == 1 ? :low : ti < N_r ? :band : :high)

	  for sgn in [1,-1]
		  # rate filtering modification (-1 for up and +1 for down)
		  if sgn > 0
			  HR = pad(HR,2length(HR))
		  else
			  HR[2:end] .= conj.(reverse(HR[2:end]))
			  HR[size(S,1)+1] = abs(HR[size(S,1)+2])
		  end

		  for (fi,scale) in enumerate(scales)
			  # scale filtering
			  HS = scale_filter(scale, size(S,1), spect_rate,
                             cm.bandonly ? :band :
                             fi == 1 ? :low : fi < size(2,1) ? :band : :high)

			  # convolve current filter with spectrogram
        Z = (st_ifft * ((HR'.*HS) .* S))
			  cr[:,(sgn==1)N_r+ti, fi, :] = view(Z,1:N_f,1:N_t)'
		  end
	  end
  end

  cr
end

function rplot(cort::CorticalRepresentation,y)
  ixs = CartesianRange(size(y))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = y[:],
                 time = times(cort,y)[at(ixs,1)][:],
                 rate = rates(cort)[at(ixs,2)][:],
                 scale = scales(cort)[at(ixs,3)][:],
                 freq_bin = at(ixs,4)[:])

  fbreaks = 2.0.^(-3:2)
  fs = freqs(cort,y)
  findices = mapslices(abs.(1000.0.*fbreaks .- fs'),2) do row
    _, i = findmin(row)
    i
  end

R"""

  library(ggplot2)

  scalestr = function(x){sprintf("Scale: %3.2f/oct",x)}
  ratestr = function(x){sprintf("Rate: %5.2f Hz",x)}
  clamp = function(x,a,b){ pmax(pmin(x,b),a) }

  df1 = $df
  df1$scale_title = factor(scalestr(df1$scale),
                           levels=scalestr($(sort(scales(cort)))))
  df1$rate_title = factor(ratestr(df1$rate),
                          levels=ratestr($(sort(rates(cort)))))
  scale = sd(df1$response)
  df1$respc = clamp(df1$response,-2*scale,2*scale)
  df1$respc[abs(df1$respc) < 2e-3] = 0.0

  ggplot(df1,aes(x=time,y=freq_bin,fill=respc)) +
    geom_raster() +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(scale_title ~ rate_title) +
    scale_fill_distiller(palette='RdBu')

"""
end
