using RCall
# without this sleep, `using RCall` can lead to a segmentation fault
# (interaction with MATLAB???)
sleep(0.25)
using Colors
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

    x = ustrip.(times(cm,data))
    y = ustrip.(freqs(cm,data))

    @series begin
      seriestype := :heatmap
      subplot := i
      (x,y,abs.(data[:,r,s,:])')
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
function askind(H,len,maxi,kind,nonorm)
  if kind in [:low,:high]
    old_sum = sum(H)
    if kind == :low
      H[1:maxi-1] = 1
    elseif kind == :high
      H[maxi+1:len] = 1
    end
    if !nonorm
      H .= H ./ sum(H) .* old_sum
    end

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

  askind(H,len,indmax(H),kind,false)
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
  HR = askind(H,len,maxi,kind,true) .* exp.(A*im)

  if rate >= 0
    HR = pad(HR,2length(HR))
	else
    HR = pad(HR,2length(HR))
		HR[2:end] .= conj.(reverse(HR[2:end]))
		HR[len+1] = abs(HR[len+2])
	end

  HR
end

(cm::CorticalModel)(s::TimedSound.Sound;usematlab=true) =
  cm(cm.aspect(s),usematlab=usematlab)
(cm::CorticalModel)(s::AbstractVector;usematlab=true) =
  cm(cm.aspect(s),usematlab=usematlab)

# TODO: use complex numbers to represent the output
# but allow plotting to show absolute value and phase???

function (cm::CorticalModel)(s_in::Matrix;usematlab=true)
  if usematlab
    s = similar(s_in,size(s_in,1),nchannels(cm.aspect))
    s[:,channels_computed(cm.aspect)] = s_in
    paras = [cm.aspect.len, cm.aspect.decay_tc,
             cm.aspect.nonlinear,cm.aspect.octave_shift,
             0.0, 0.0, cm.bandonly ? 1.0 : 0.0]
    if !all(indexin(-abs.(cm.rates),cm.rates) .> 0)
      error("Missing negative rates. Cannot use matlab implementation."*
            " Set usematlab=false.")
    end
    if any(isnan.(cm.rates)) || any(isnan.(cm.scales))
      error("NaN rate and/or scale not supproted by matlab implementation."*
            " Set usematlab=false.")
    end
    y = mat"aud2cor($s,$paras,$(unique(abs.(cm.rates))),$(cm.scales),'tmpxxx',0)"

    y = permutedims(y,[3,2,1,4])
    rs = Main.rates(cort)
    order = sortperm([.-reverse(rs[rs .> 0]); reverse(rs[rs .> 0])])

    return y[:,order,:,channels_computed(cm.aspect)]
  end

  s = s_in
  warn("Julia cortical result differs from NSL toolbox result!!!")

  rates = cm.rates
  scales = cm.scales
  N_t, N_f = map(n -> nextprod([2,3,5],n),size(s)) # TODO: change to [2,3,5]
  N_r, N_s = length(rates), length(scales)

  # spatial-frequency represention of the spectrogram (i.e. 2D fft of s).
  S1 = fft(pad(s,(N_t,2N_f)),2)
  S = fft(pad(S1[:,1:N_f],(2N_t,N_f)),1)
  yt = s

  t_ifft = plan_ifft(S,1)
  f_ifft = plan_ifft(Array{eltype(s)}(size(s,1),2N_f),2)

  S1 = S1[:,1:N_f]
  s1_ifft = plan_ifft(S1,1)

  cr = zeros(Complex{eltype(s)}, size(s,1), N_r, N_s, size(s,2))
  rmin,rmax = extrema(rates)
  smin,smax = extrema(scales)
  @showprogress "Cortical Simulation: " for (ri,rate) in enumerate(rates)
    z_t = if isnan(rate)
      # do not filter by rate
      (t_ifft * S)[1:size(s,1),:]
    else
	    HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                       cm.bandonly ? :band :
                       rate == rmin ? :low : rate < rmax ? :band : :high)

      # apply the rate filter
      (t_ifft * (HR.*S))[1:size(s,1),:]
    end

    for (si,scale) in enumerate(scales)
      if isnan(scale)
        # do not filter by scale
        z = t_ifft * (HR .* S1)
        cr[:, ri, si, :] = view(z,indices(s)...)
      else
			  HS = scale_filter(scale, N_f, spect_rate,
                          cm.bandonly ? :band :
                          scale == smin ? :low : scale < smax ? :band : :high)

			  # apply the scale filter
        z = f_ifft*(pad(z_t.*HS',size(z_t,1),2N_f))

			  cr[:, ri, si, :] = view(z,indices(s)...)
      end
		end
  end

  cr
end

function Base.inv(cm::CorticalModel,cr_in::Array{T,4};
                  usematlab=true,norm=0.9) where T
  if usematlab
    cr = similar(cr_in,(size(cr_in,1,2,3)...,nchannels(cm.aspect)))
    cr[:,:,:,channels_computed(cm.aspect)] = cr_in

    if !all(indexin(-abs.(cm.rates),cm.rates) .> 0)
      error("Missing negative rates. Cannot use matlab implementation."*
            " Set usematlab=false.")
    end
    if any(isnan.(cm.rates)) || any(isnan.(cm.scales))
      error("NaN rate and/or scale not supproted by matlab implementation."*
            " Set usematlab=false.")
    end
    m_cr = permutedims(cr,[3,2,1,4])

    # re-order by rates
    rs = Main.rates(cort)
    order = sortperm([.-reverse(rs[rs .> 0]); reverse(rs[rs .> 0])])
    inv_order = sortperm(order)
    m_cr = m_cr[:,inv_order,:,:]

    cm(cr[:,1,1,channels_computed(cm.aspect)],usematlab=true) # this ensures 'tmpxxx' is right, hacky but it works
    mat"s = cor2aud('tmpxxx',$m_cr);"
    s = mat"s + 0"
    s = 2.*max.(real.(s),0)
    return s[:,channels_computed(cm.aspect)]
  else
    cr = cr_in

    warn("Julia implementation of inverse is very poor!")
    # TODO: look at intermediate outcomes to troubleshoot the implementation

    rates = cm.rates
    scales = cm.scales
    N_t, N_f = map(n -> nextprod([2,3,5],n),size(cr,1,4)) # TODO: change to [2,3,5]
    N_r, N_s = size(cr,2,3)

    z = zeros(T,2N_t,2N_f)
    z_cum = zeros(z)
    h_cum = zeros(real(T),size(z)...)
    st_fft = plan_fft(z_cum)

    rmin,rmax = extrema(rates)
    smin,smax = extrema(scales)
    @showprogress "Inverting Cortical Simulation: " for (ri,rate) in enumerate(rates)
      if isnan(rate) continue end
	    # rate filtering
      HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                       cm.bandonly ? :band :
                       rate == rmin ? :low : rate < rmax ? :band : :high)

      for (si,scale) in enumerate(scales)
        if isnan(rate) continue end
			  # scale filtering
			  HS = scale_filter(scale, N_f, spect_rate,
                          cm.bandonly ? :band :
                          scale == smin ? :low : scale < smax ? :band : :high)

        z[indices(cr,1),indices(cr,4)] = cr[:,ri,si,:]

        Z = st_fft * z
        h = HR.*[HS; zeros(HS)]'
        h_cum .+= abs2.(h)
        z_cum .+= h .* Z
      end
    end

    h_cum[:,1] .*= 2
    old_sum = sum(h_cum)
    h_cum .= norm.*h_cum + (1 .- norm).*maximum(h_cum)
    h_cum .*= old_sum ./ sum(h_cum)
    z_cum ./= h_cum

    s = (st_fft \ z_cum)[indices(cr,1),indices(cr,4)]
    max.(real.(2.*s),0)
  end
end

rplot(cort::CorticalModel,y::AbstractVector) = rplot(cort,cort(y))

# function number_to_color(x::Array{<:Complex})
#   phase_map = cmap("C9")
#   norm(x) = x ./ maximum(x)
#   tocolor(x) = phase_map[floor(Int,x*(length(phase_map)-1)+1)]

#   phase_col = Lab.(tocolor.((angle.(x) .+ π)/2π))
#   RGB.(weighted_color_mean.(norm(log.(1 .+ abs.(x))),
#                             phase_col,Lab(colorant"lightgray")))
# end

# function number_to_color(x::Array{<:Real})
#   colormap = cmap("L3")
#   tocolor(x) = colormap[floor(Int,x*(length(colormap)-1)+1)]

#   RGB.(tocolor.((x .- minimum(x)) ./ (maximum(x) - minimum(x))))
# end

function nearin(xin,xs)
  _,inds = findmin(abs.(vec(xin) .- vec(xs)'),2)
  cols = map(ii -> ii[2],CartesianRange((length(xin),length(xs))))
  cols[inds[:,1]]
end

function rplot(cort::CorticalModel,y;rates=cort.rates,scales=cort.scales)
  rindices = nearin(rates,cort.rates)
  sindices = nearin(scales,cort.scales)

  if rates != cort.rates || scales != cort.scales
    @show rindices
    @show sindices
  end

  y = y[:,rindices,sindices,:]

  ixs = CartesianRange(size(y))
  at(ixs,i) = map(x -> x[i],ixs)

  colormap = "#".*hex.(RGB.(cmap("C6")))

  df = DataFrame(r_phase = angle.(vec(y)),
                 r_amp = abs.(vec(y)),
                 time = ustrip(vec(times(cort,y)[at(ixs,1)])),
                 rate = vec(rates[at(ixs,2)]),
                 scale = vec(scales[at(ixs,3)]),
                 freq_bin = vec(at(ixs,4)))

  fbreaks,findices = freq_ticks(cort.aspect,y[:,1,1,:])

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

  ggplot(df1,aes(x=time,y=freq_bin,fill=r_phase,alpha=r_amp)) +
    geom_raster() +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(scale_title ~ rate_title) +
    scale_fill_gradientn(colors=$colormap,limits=c(-pi-0.1,pi+0.01),
                         breaks=c(-pi,0,pi),
                         labels=c(expression(-pi),expression(0),
                                  expression(+pi)),
                                  name = expression(phi)) +
    scale_alpha_continuous(range=c(0,1),name="Amplitude")

"""
end

function plot_scales(cort,m,range=nothing)
  sm = m[:,1,:,1]
  iis = collect(CartesianRange(size(sm)))
  df = DataFrame(level = vec(sm),
                 time = vec(map(x -> ustrip(times(cort,m)[x[1]]),iis)),
                 scale = vec(map(x -> ustrip(scales(cort)[x[2]]),iis)))

  @show unique(df[:scale])
R"""
  library(RColorBrewer)
  pal = brewer.pal(5,'Reds')[2:5]
  p = ggplot($df,aes(x=time,y=level,group=scale,
                 color=factor(round(scale,1)),linetype=factor(round(scale,1)))) +
    geom_line(color='black',linetype='solid',size=1.2) + geom_line() +
    scale_color_manual(values = rep(pal,each=3)[1:11],name="Scale") +
    scale_linetype_manual(values =
      rep(c("dotdash","longdash","solid"),4)[1:11],name="Scale")
"""
  if range != nothing
R"""
      p = p + coord_cartesian(ylim=c($(first(range)),$(last(range))))
"""
  end

  R"p"
end

function plot_scales2(cort,data,range=nothing)
  data = data[:,1,:,1]
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(cort,data)[at(ixs,1)])),
                 scale_bin = vec(at(ixs,2)))

  sbreaks = 1:2:length(scales(cort))
  slabs = string.(round(scales(cort)[sbreaks],1)).*" cyc/oct"

  @show sc
  @show unique(df[:scale_bin])
R"""

  library(ggplot2)

  print($sc)

  ggplot($df,aes(x=time,y=scale_bin,fill=response)) +
    geom_raster() +
    scale_y_continuous(labels=$slabs,breaks=$sbreaks) +
    ylab('Frequency (Hz)') + xlab('Time (s)') +
    scale_fill_distiller(palette='Reds',direction=1)

"""
end
