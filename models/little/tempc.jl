using Parameters
using ProgressMeter
include("cortical.jl")
include("online_pca.jl")

struct TCAnalysis
  upstream::CorticalModel
  ncomponents::Int
  rate::Seconds{Float64}
  method::Symbol
end

Δt(tc::TCAnalysis) = Δt(tc.upstream)
times(tc::TCAnalysis,x) = times(tc.upstream,x)

(tc::TCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
(tc::TCAnalysis)(x::AbstractMatrix) = tc(tc.upstream(x))

TCAnalysis(upstream,ncomponents,rate=1s;method=:pca) =
  TCAnalysis(upstream,ncomponents,rate,method)

# alternative: I could have a different
# set of eigenseries for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (tc::TCAnalysis)(x)
  if tc.method == :pca
    C = EigenSeries(eltype(x),size(x,1),prod(size(x,3,4)),tc.ncomponents,Δt(tc))
    window_len = round(Int,tc.rate/Δt(tc))
    @showprogress "Temporal Coherence Analysis: " for t in indices(x,1)
      x_ts = x[max(1,t-window_len):t,:,:,:]
      x_t = reshape(x_ts,prod(size(x_ts,1,2)),:)

      # TODO: change to approximate: C_t = (1-dt)C_t + x*x'*dt ??
      # (by weighting older samples)
      n = min(size(x_t,1),tc.ncomponents)
      sv, = svds(x_t,nsv=n)
      λ = zeros(eltype(x),tc.ncomponents)
      u = zeros(eltype(x),size(x_t,2),tc.ncomponents)

      λ[1:n] = sv[:S].^2 ./ size(x_t,1)
      u[:,1:n] = sv[:V]
      var = mean(abs2.(x_t),1)
      C[t] = EigenSpace(sv[:V],(sv[:S]).^2 / size(x_t,1),var)
    end

    C
  else
    error("No method named $(tc.method)")
  end
end

fusion_signal(tc::TCAnalysis,C::EigenSeries) = fusion_signal(tc.upstream(x),C)
function fusion_signal(tc::TCAnalysis,C::EigenSeries)
  vec(first.(eigvals.(C)) ./ sum.(var.(C)))
end

mask(tc::TCAnalysis,C::EigenSpace,x,phase;component=1) =
    mask(tc,C,tc.upstream(x),phase,component=component)
function mask(tc::TCAnalysis,C::EigenSpace,x::Array{T,4},
              phase;component=1) where T
  m = real.(eigvecs(C)[:,component] .* exp(phase*im))
  mr = reshape(m,size(x,3,4)...)
  mr ./= maximum(abs.(mr))
  mr .= max.(0,mr)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = mr.*x[ii,:,:]
  end
  y
end


rplot(tc::TCAnalysis,x::TimedSound.Sound;kwds...) = rplot(tc,tc(x);kwds...)

function rplot(tc::TCAnalysis,λ::Vector)
  @assert all(imag.(λ) .== 0) "Can't plot complex eigenvalues"
  λ = real.(λ)
  df = DataFrame(value = sort(λ,rev=true),index = collect(eachindex(λ)))

R"""
  library(ggplot2)

  ggplot($df,aes(x=index,y=value)) + geom_bar(stat='identity')
"""
end

# TODO: improve this plot by making it relative to variance (ala fusion signal)
# and then make plot_resps a version of this using an array of eigenseries
function rplot(tc::TCAnalysis,C::EigenSeries;n=ncomponents(C))
  λ = eigvals(C)
  λ = λ[:,sortperm(abs.(λ[max(1,end-10),:]),rev=true)]
  ii = CartesianRange(size(λ))
  at(i) = map(ii -> ii[i],ii)

  prop_complex = mean(imag.(λ))
  if prop_complex > 0
    warn("$(round(100prop_complex,1))% of eigenvalues were "*
         "complex (using absolute value).")
  end

  df = DataFrame(value = vec(abs.(λ) ./ sum.(var.(C))),
                 time = vec(ustrip(at(1) * Δt(spect))),
                 component = vec(at(2)))

  df = df[df[:component] .<= n,:]

R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=value,color=factor(component),group=component)) +
    geom_line() + scale_color_brewer(palette='Set1',name='Component') +
    xlab('Time (s)') + ylab('Value')
"""
end

function rplot(tc::TCAnalysis,C::EigenSpace;n=ncomponents(C),showvar=true)
  λ = abs.(eigvals(C))
  order = sortperm(λ,rev=true)
  λ = λ[order]
  u = eigvecs(C)[:,order]
  u = u[:,1:min(n,end)]
  if showvar; u = [u var(C)]; end
  u = reshape(u,length(scales(tc.upstream)),:,size(u,2))
  ii = CartesianRange(size(u))
  at(i) = vec(map(ii -> ii[i],ii))
  function title(n)
    if n <= length(λ)
      # nstr = @sprintf("%02d",n)
      "lambda[$(n)] == $(round(λ[n],1))"
      # "Lmb_$nstr = $(round(λ[n],3))"
    else
      "sigma^2 == $(round(sum(var(C)),1))"
      # "Variance = $(sum(var(C)))"
    end
  end

  colormap = "#".*hex.(RGB.(cmap("C6")))

  df = DataFrame(r_phase = angle.(vec(u)) ,
                 r_amp = abs.(vec(u)),
                 scale_index = at(1),
                 freq_bin = at(2),
                 component = at(3),
                 component_title = title.(at(3)))

  sindices = 1:2:length(scales(tc.upstream))
  sbreaks = scales(tc.upstream)[sindices]
  fbreaks,findices = freq_ticks(tc.upstream.aspect,u[:,:,1])

R"""

  library(ggplot2)

  ggplot($df,aes(x=scale_index,y=freq_bin,fill=r_phase,alpha=r_amp)) +
    geom_raster() + facet_wrap(~component_title,labeller=label_parsed) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    scale_x_continuous(breaks=$sindices,labels=$sbreaks) +
    ylab('Frequency (Hz)') + xlab('Scale (cycles/octave)') +

    scale_fill_gradientn(colors=$colormap,limits=c(-pi-0.01,pi+0.01),
                         breaks=c(-pi,0,pi),
                         labels=c(expression(-pi),expression(0),
                                  expression(+pi)),
                         name = expression(phi))+
    scale_alpha_continuous(range=c(0,1),name="Amplitude")

"""

end

function rplot(tempc::TCAnalysis,C::Array{<:EigenSeries},
               label=:index => 1:length(C))
  x = [fusion_signal(tempc,Ci) for Ci in C]
  df = DataFrame(resp = vcat((real.(xi) for xi in x)...),
                 var = vcat((fill(label[2][i],length(x[i]))
                             for i in eachindex(x))...),
                 time = vcat((ustrip.(eachindex(xi) * Δt(spect))
                              for xi in x)...))

R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=resp,group=factor(var),color=factor(var))) +
    geom_line() +
    scale_color_brewer(palette='Set1',name=$(string(label[1]))) +
    coord_cartesian(ylim=c(1.1,0)) + ylab('lambda / var(x)') +
    xlab('time (s)')
"""
end
