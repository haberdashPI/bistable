using RCall
# without this sleep, `using RCall` can lead to a segmentation fault
# (interaction with MATLAB???)
sleep(0.25)
using Colors
using PerceptualColourMaps
using RecipesBase
using DataFrames

@static if USING_MATLAB using MATLAB end

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

const loadloaded = fill(false)
function CorticalModel(aspect::AuditorySpectrogram;
                       rates=sort([-2.^(1:0.5:5); 2.^(1:0.5:5)]),
                       scales=2.^(-2:0.5:3),
                       bandonly=false)
  CorticalModel(aspect,rates,scales,bandonly)
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
  if kind == :band
    H
  else
    old_sum = sum(H)
    if kind == :low
      H[1:maxi-1] = 1
    elseif kind == :high
      H[maxi+1:len] = 1
    else
      error("Unexpected filter kind '$kind'.")
    end
    if !nonorm
      H .= H ./ sum(H) .* old_sum
    end

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
function rate_filter(rate,len,spect_len,kind,use_conj=false)
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

  if use_conj
    HR = conj.(HR)
  end

  if rate >= 0
    HR = pad(HR,2length(HR))
	else
    HR = pad(HR,2length(HR))
		HR[2:end] .= conj.(reverse(HR[2:end]))
		HR[len+1] = abs(HR[len+2])
	end

  HR
end

(cm::CorticalModel)(s::Sound;usematlab=false) =
  cm(cm.aspect(s),usematlab=usematlab)
(cm::CorticalModel)(s::AbstractVector;usematlab=false) =
  cm(cm.aspect(s),usematlab=usematlab)

# TODO: use complex numbers to represent the output
# but allow plotting to show absolute value and phase???

function (cm::CorticalModel)(s_in::Matrix;usematlab=false,progressbar=true)
  if usematlab
    @static if USING_MATLAB
      if !loadloaded[]
        loadloaded[] = true
        mat"loadload;"
      end
      s = similar(s_in,(size(s_in,1),nchannels(cm.aspect)))
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
      orates = sort(unique(abs.(cm.rates)))
      oscales = sort(unique(cm.scales))
      y = mat"aud2cor($s,$paras,$orates,$oscales,'tmpxxx',0)"

      y = permutedims(y,[3,2,1,4])
      rs = Main.rates(cort)
      order = sortperm([.-rs[rs .> 0]; rs[rs .> 0]])

      return y[:,order,:,channels_computed(cm.aspect)]
    else
      error("Matlab implementation not supported. Set `USING_MATLAB = true`.")
    end
  end

  s = s_in

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
  rmin,rmax = extrema(abs.(rates))
  smin,smax = extrema(scales)
  progress = Progress(length(rates)*length(scales),
                      desc="Cortical Simulation: ",
                      dt = progressbar ? 1.0 : Inf)
  for (ri,rate) in enumerate(rates)
    z_t = if isnan(rate)
      # do not filter by rate
      (t_ifft * S)[1:size(s,1),:]
    else
      HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                       cm.bandonly ? :band :
                       abs(rate) == rmin ? :low :
                       abs(rate) < rmax ? :band : :high)

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

      next!(progress)
		end
  end

  cr
end

function Base.inv(cm::CorticalModel,cr_in::Array{T,4};
                  usematlab=false,norm=0.9,progressbar=true) where T
  if usematlab
    @static if USING_MATLAB
      if !loadloaded[]
        loadloaded[] = true
        mat"loadload;"
      end

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
      order = sortperm([.-rs[rs .> 0]; rs[rs .> 0]])
      inv_order = sortperm(order)
      m_cr = m_cr[:,inv_order,:,:]

      cm(cr[:,1,1,channels_computed(cm.aspect)],usematlab=true) # this ensures 'tmpxxx' is right, hacky but it works
      mat"s = cor2aud('tmpxxx',$m_cr);"
      s = mat"s + 0"
      s = 2.*max.(real.(s),0)
      return s[:,channels_computed(cm.aspect)]
    else
      error("Matlab implementation not supported. Set `USING_MATLAB = true`.")
    end
  else
    cr = cr_in

    rates = cm.rates
    scales = cm.scales
    N_t, N_f = map(n -> nextprod([2],n),size(cr,1,4)) # TODO: change to [2,3,5]
    N_r, N_s = size(cr,2,3)

    z = zeros(T,2N_t,2N_f)
    z_cum = zeros(z)
    h_cum = zeros(real(T),size(z)...)
    st_fft = plan_fft(z_cum)

    rmin,rmax = extrema(abs.(rates))
    smin,smax = extrema(scales)
    dt = progressbar ? 1.0 : Inf
    progress = Progress(length(rates)*length(scales),
                        desc="Cortical Inversion: ",
                        dt = progressbar ? 1.0 : Inf)
    for (ri,rate) in enumerate(rates)
      if isnan(rate)
        next!(progresse)
        continue
      end
	    # rate filtering
      HR = rate_filter(rate, N_t, 1000 / cm.aspect.len,
                       cm.bandonly ? :band :
                       abs(rate) == rmin ? :low :
                       abs(rate) < rmax ? :band : :high,
                       true)

      for (si,scale) in enumerate(scales)
        if isnan(rate)
          next!(progress)
          continue
        end
			  # scale filtering
			  HS = scale_filter(scale, N_f, spect_rate,
                          cm.bandonly ? :band :
                          scale == smin ? :low : scale < smax ? :band : :high)

        z[indices(cr,1),indices(cr,4)] = cr[:,ri,si,:]

        Z = st_fft * z
        h = HR.*[HS; zeros(HS)]'
        h_cum .+= abs2.(h)
        z_cum .+= h .* Z

        next!(progress)
      end
    end

    h_cum[:,1] .*= 2
    old_sum = sum(h_cum[:,indices(cr,4)])
    h_cum .= norm.*h_cum + (1 .- norm).*maximum(h_cum)
    h_cum .*= old_sum ./ sum(h_cum[:,indices(cr,4)])
    z_cum ./= h_cum

    s = (st_fft \ z_cum)[indices(cr,1),indices(cr,4)]
    max.(real.(2.*s),0)
  end
end
