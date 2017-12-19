# TODO: figure out if the segfault comes from the combination of R and MATLAB
using Colors
using RCall
using MATLAB
using PerceptualColourMaps
using RecipesBase
using DataFrames
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
  mat"loadload;"
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
    ones(len)
  end

  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)

  askind(H,len,indmax(H),kind)
end

# create the temporal-rate filter (filter along temporal axis)
function rate_filter(rate,len,spect_len,kind)
  if isnan(rate)
    return ones(2len)
  end

  t = (0:len-1)./spect_len .* abs(rate)
  h = @. abs(rate) * sin(2π*t) * t^2 * exp(-3.5t)
  h .-= mean(h)

  H0 = view(fft(pad(h,2len)),1:len)
  A = angle.(H0)
  H = abs.(H0)

  maxH,maxi = findmax(H)
  H ./= maxH
  HR = askind(H,len,maxi,kind) .* exp.(A*im)

  if rate >= 0
    HR = pad(HR,2length(HR))
	else
    HR = pad(HR,2length(HR))
		HR[2:end] .= conj.(reverse(HR[2:end]))
		HR[len+1] = abs(HR[len+2])
	end

  HR
end

(cm::CorticalModel)(s::TimedSound.Sound,usematlab=true) =
  cm(cm.aspect(s),usematlab)
(cm::CorticalModel)(s::AbstractVector,usematlab=true) =
  cm(cm.aspect(s),usematlab)

# TODO: use complex numbers to represent the output
# but allow plotting to show absolute value and phase???

function (cm::CorticalModel)(s::Matrix,usematlab=true)
  if usematlab
    paras = [cm.aspect.len, cm.aspect.decay_tc,
             cm.aspect.nonlinear,cm.aspect.octave_shift]
    if !all(indexin(-abs.(cm.rates),cm.rates) .> 0)
      error("Missing negative rates. Cannot use matlab implementation.")
    end
    y = mat"aud2cor($s,$paras,$(unique(abs.(cm.rates))),$(cm.scales),'tmpxxx',0)"
    return permutedims(y,[3,2,1,4])
  end

  warn("Julia cortical implementation is not fully functional!!!")

  rates = cm.rates
  scales = cm.scales
  N_t, N_f = map(n -> nextprod([2],n),size(s)) # TODO: change to [2,3,5]
  N_r, N_s = length(rates), length(scales)

  # spatial-frequency represention of the spectrogram (i.e. 2D fft of s).
  S1 = fft(pad(s,(N_t,2N_f)),2)
  S = fft(pad(S1[:,1:N_f],(2N_t,N_f)),1)

  t_ifft = plan_ifft(S,1)
  f_ifft = plan_ifft(S1,2)

  S1 = S1[:,1:N_f]
  s1_ifft = plan_ifft(S1,1)

  cr = zeros(Complex{eltype(s)}, size(s,1), N_r, N_s, size(s,2))
  rmin,rmax = extrema(rates)
  smin,smax = extrema(scales)
  @showprogress "Cortical Simulation " for (ri,rate) in enumerate(rates)
	  # rate filtering
	  HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                     cm.bandonly ? :band :
                     rate == rmin ? :low : rate < rmax ? :band : :high)

    z_t = (t_ifft * (HR.*S))[1:N_t,:]

    for (si,scale) in enumerate(scales)
      if isnan(scale)
        z = t_ifft * (HR .* S1)
        cr[:, ri, si, :] = view(z,indices(s)...)
      else
			  # scale filtering
			  HS = scale_filter(scale, N_f, spect_rate,
                          cm.bandonly ? :band :
                          scale == smin ? :low : scale < smax ? :band : :high)

			  # convolve current filter with spectrogram
        z = f_ifft*(pad(HS'.*z_t,N_t,2N_f))
			  cr[:, ri, si, :] = view(z,indices(s)...)
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

function number_to_color(x::Array{<:Complex})
  phase_map = cmap("C9")
  norm(x) = x ./ maximum(x)
  tocolor(x) = phase_map[floor(Int,x*(length(phase_map)-1)+1)]

  phase_col = Lab.(tocolor.((angle.(x) .+ π)/2π))
  RGB.(weighted_color_mean.(norm(log.(1 .+ abs.(x))),
                            phase_col,Lab(colorant"lightgray")))
end

function number_to_color(x::Array{<:Real})
  colormap = cmap("L3")
  tocolor(x) = colormap[floor(Int,x*(length(colormap)-1)+1)]

  RGB.(tocolor.((x .- minimum(x)) ./ (maximum(x) - minimum(x))))
end

function rplot(cort::CorticalModel,y;rates=cort.rates,scales=cort.scales)
  rindices = indexin(rates,cort.rates)
  sindices = indexin(scales,cort.scales)

  if rates != cort.rates || scales != cort.scales
    @show rindices
    @show sindices
  end

  y = y[:,rindices,sindices,:]

  ixs = CartesianRange(size(y))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(r_color = "#".*hex.(number_to_color(y[:])),
                 time = times(cort,y)[at(ixs,1)][:],
                 rate = rates[at(ixs,2)][:],
                 scale = scales[at(ixs,3)][:],
                 freq_bin = at(ixs,4)[:])

  @show head(df)

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
