using Colors
using PerceptualColourMaps
using RecipesBase
using DataFrames
using RCall
include("audio_spect.jl")

struct CorticalModel
  aspect::AuditorySpectrogram
  rates::Vector{Float64}
  scales::Vector{Float64}
  bandonly::Bool
end

scales(cm::CorticalModel) = cm.scales
rates(cm::CorticalModel) = cm.rates

freqs(cm::CorticalModel,data::Array{T,4}) where T =
  @views freqs(cm.aspect,data[:,1,1,:])
freqs(cm::CorticalModel,data::Array{T,2}) where T =
  @views freqs(cm.aspect,data)

times(cm::CorticalModel,data::Array{T,4}) where T =
  @views times(cm.aspect,data[:,1,1,:])
times(cm::CorticalModel,data::Array{T,2}) where T =
  @views times(cm.aspect,data)

function CorticalModel(aspect::AuditorySpectrogram;
                       rates=sort([-2.^(1:0.5:5); 2.^(1:0.5:5)]),
                       scales=2.^(-2:0.5:3),
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
  if isnan(scale)
    ones(2len)
  end

  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)
  _,maxi = findmax(H)

  H = askind(H,len,maxi,kind)
  H = [zeros(H); H]
end

# create the temporal-rate filter (filter along temporal axis)
function rate_filter(rate,len,spect_len,kind)
  if isnan(rate)
    return ones(2len)
  end

  t = (0:len-1)./spect_len .* abs(rate)
  h = @. abs(rate) * t^3 * exp(-4t) * cos(2π*t)
  h .-= mean(h)

  H0 = view(fft(pad(h,2len)),1:len)
  A = angle.(H0)
  H = abs.(H0)

  _, maxi = findmax(H)
  H ./= H[maxi]
  HR = askind(H,len,maxi,kind) .* exp.(A*im)

  if rate >= 0
    HR = pad(HR,2length(HR))
	else
    HR = pad(HR,2length(HR))
    HR[1] = 0
		HR[2:end] .= conj.(reverse(HR[2:end]))
		HR[len+1] = abs(HR[len+2])
	end

  HR
end

(cm::CorticalModel)(s::AbstractVector) = cm(cm.aspect(s))

# TODO: use complex numbers to represent the output
# but allow plotting to show absolute value and phase???

function (cm::CorticalModel)(s::Matrix)
  nexteven(x) = iseven(x) ? x : x+1

  rates = cm.rates
  scales = cm.scales
  N_t, N_f = nexteven.(size(s))
  N_r, N_s = length(rates), length(scales)

  # spatial-frequency represention of the spectrogram (i.e. 2D fft of s).

  # implementation note: we transpose s because rfft halves the first dimension
  # and we want frequency not time to be the halved dimension (since the rate
  # filters, over time, need both negative and positive direction and that
  # invovles manipulating negative frequencies)
  S = fft(pad(s,nexteven.(size(s))))
  S1 = fft(pad(s,nexteven(size(s,1)),size(s,2)),1)
  st_ifft = plan_ifft(S)
  t_ifft = plan_ifft(S1,1)

  cr = zeros(Complex{eltype(s)}, N_t, N_r, N_s, N_f)
  rmin,rmax = extrema(rates)
  smin,smax = extrema(scales)
  for (ri,rate) in enumerate(rates)
	  # rate filtering
	  HR = rate_filter(rate, N_t>>1, 1000 / cm.aspect.len,
                     cm.bandonly ? :band :
                     rate == rmin ? :low : rate < rmax ? :band : :high)

    for (si,scale) in enumerate(scales)
      if isnan(scale)
        Z = t_ifft * (HR .* S1)
        cr[:, ri, si, :] = view(Z,indices(s)...)
      else
			  # scale filtering
			  HS = scale_filter(scale, N_f>>1, spect_rate,
                          cm.bandonly ? :band :
                          scale == smin ? :low : scale < smax ? :band : :high)

			  # convolve current filter with spectrogram
        Z = st_ifft*(HR.*HS' .* S)
			  cr[:, ri, si, :] = view(Z,indices(s)...)
      end
		end
  end

  cr
end


#========================================
# TODO: code review/debug/test
struct InvCorticalModel
  cm::CorticalModel
end
Base.inv(cm::CorticalModel) = InvCorticalModel(cm)

function (cmi::InvCorticalModel)(cr::Array{T,4}) where T
  cm = cmi.cm
  rates = cm.rates
  scales = cm.scale
  N_t, N_f = size(s)
  N_r, N_s = length(rates), length(scales)

  z = zeros(T,2N_t,2N_f)
  z_cum = zeros(z) # ????
  h_cum = zeros(z)
  st_fft = plan_fft(spect)

  for (ri,rate) in enumerate(rates)
    if isnan(rate) break end
	  # rate filtering
	  HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                     cm.bandonly ? :band :
                     rate == rmin ? :low : rate < rmax ? :band : :high)

    for (si,scale) in enumerate(scales)
      if isnan(rate) break end
			# scale filtering
			HS = scale_filter(scale, size(S,1), spect_rate,
                        cm.bandonly ? :band :
                        scale == smin ? :low : scale < smax ? :band : :high)

      z[1:N_t,1:N_f] = cr[:,ri,si,:]

      Z = st_fft * z
      h = HR.*HS'
      h_cum .+= h .* conj.(h)
      z_cum .+= h .* Z
    end
  end

  h_cum[:,1] .*= 2
  old_sum = sum(h_cum)
  h_cum .= cmi.norm.*h_cum + (1 .- cmi.norm).*maximum(h_cum)
  h_cum .*= old_sum ./ sum(h_cum)

  z_cum[1:2Nt,1:N_f] ./= h_cum[1:2n_t,1:N_f]

  2(st_fft \ z_cum)
end
========================================#

rplot(cort::CorticalModel,y::AbstractVector) = rplot(cort,cort(y))

function rplot(cort::CorticalModel,y;rates=cort.rates,scales=cort.scales)
  rindices = indexin(rates,cort.rates)
  sindices = indexin(scales,cort.scales)

  @show rindices
  @show sindices

  y = y[:,rindices,sindices,:]

  ixs = CartesianRange(size(y))
  at(ixs,i) = map(x -> x[i],ixs)
  norm(x) = x ./ maximum(x)

  phase_map = cmap("C9")
  tocolor(x) = phase_map[floor(Int,x*(length(phase_map)-1)+1)]
  phase_col = Lab.(tocolor.((angle.(y[:]) .+ π)/2π))
  col = weighted_color_mean.(norm(log.(1 + abs.(y[:]))),
                             phase_col,Lab(colorant"lightgray"))

  df = DataFrame(r_color = "#".*hex.(col),
                 time = times(cort,y)[at(ixs,1)][:],
                 rate = rates[at(ixs,2)][:],
                 scale = scales[at(ixs,3)][:],
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
                           levels=scalestr($(sort(scales))))
  df1$rate_title = factor(ratestr(df1$rate),
                          levels=ratestr($(sort(rates))))

  ggplot(df1,aes(x=time,y=freq_bin,fill=r_color)) +
    geom_raster() +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(scale_title ~ rate_title) +
    scale_fill_identity() # scale_fill_distiller(palette='RdBu')

"""
end
