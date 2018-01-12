using Parameters
using ProgressMeter
include("cortical.jl")
include("online_pca.jl")

struct TCAnalysis
  cort::CorticalModel
  ncomponents::Int
  window::Seconds{Float64}
  method::Symbol
  frame_len::Int
end

Δt(tc::TCAnalysis) = Δt(tc.cort)*tc.frame_len
times(tc::TCAnalysis,x) = times(tc.cort,x)[1:tc.frame_len:end]

(tc::TCAnalysis)(x::AbstractVector) = tc(tc.cort(x))
(tc::TCAnalysis)(x::AbstractMatrix) = tc(tc.cort(x))

TCAnalysis(cort,ncomponents;window=1s,method=:pca,frame_len=10ms) =
  TCAnalysis(cort,ncomponents,window,method,
             max(1,floor(Int,frame_len/Δt(cort))))

windowing(x,dim;length=nothing,step=nothing) =
  (max(1,t-length):t for t in indices(x,dim)[1:step:end])

windowlen(tc::TCAnalysis) = round(Int,tc.window/Δt(tc))
nunits(tc::TCAnalysis,x) = prod(size(x,3,4))
ncomponents(tc::TCAnalysis) = tc.ncomponents

# alternative: I could have a different
# set of eigenseries for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (tc::TCAnalysis)(x)
  if tc.method == :pca
    windows = enumerate(windowing(x,1;length=windowlen(tc),step=tc.frame_len))
    C = EigenSeries(eltype(x),length(windows),nunits(tc,x),
                    ncomponents(tc),Δt(tc))

    @showprogress "Temporal Coherence Analysis: " for (i,w_inds) in windows
      window = x[w_inds,:,:,:]
      x_t = reshape(window,prod(size(window,1,2)),:)

      n = min(size(x_t,1),ncomponents(tc))
      sv, = svds(x_t,nsv=n)

      λ = zeros(eltype(x),ncomponents(tc))
      u = zeros(eltype(x),size(x_t,2),ncomponents(tc))

      λ[1:n] = sv[:S].^2 ./ size(x_t,1)

      u[:,1:n] = sv[:V]
      var = mean(abs2.(x_t),1)
      C[i] = EigenSpace(sv[:V],(sv[:S]).^2 / size(x_t,1),var)
    end

    C
  else
    error("No method named $(tc.method)")
  end
end

fusion_signal(tc::TCAnalysis,C::EigenSeries) = fusion_signal(tc.cort(x),C)
function fusion_signal(tc::TCAnalysis,C::EigenSeries)
  vec(first.(eigvals.(C)) ./ sum.(var.(C)))
end

const phase_resolution = 128

mask(tc::TCAnalysis,C,x;kwds...) = mask(tc,C,tc.cort(x);kwds...)
function mask(tc::TCAnalysis,C::EigenSpace,x::Array{T,4};
              phase=nothing,component=1) where T

  M(phase) = max.(0,real.(eigvecs(C)[:,component] .* exp.(phase*im)))

  # if there is no phase specified, choose the maximum energy mask
  phase = if phase !== nothing; phase; else
    phases = linspace(-π,π,phase_resolution+1)[1:end-1]
    phases[indmax(sum.(M.(phases)))]
  end

  m = reshape(M(phase),size(x,3,4)...)
  m ./= maximum(m)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = m.*x[ii,:,:]
  end
  y
end

function mask(tc::TCAnalysis,C::EigenSeries,x::Array{T,4}) where T
  y = zero(x)
  norm = fill(0.0,size(y,1))

  windows = enumerate(windowing(x,1;length=windowlen(tc),step=tc.frame_len))
  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask(tc,C[i],window)
    y[w_inds,:,:,:] .+= m
    norm[w_inds] += 1
  end

  y ./= max.(1e-10,norm)
end

# in general, we have to solve the problem by
# computing for each time point, over the same
# window used to compute the component.

# NOTE: that doesn't seem to really work all that well

# a potentially simple hueristic for ABA stimuli:
# compute energy at the maximum energy phase
# and compare to the energy at max(phase) + π

# make a plot of both these estimates
# to see how they compare

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
                 time = vec(ustrip(at(1) * Δt(tc))),
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
  u = reshape(u,length(scales(tc.cort)),:,size(u,2))
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

  sindices = 1:2:length(scales(tc.cort))
  sbreaks = scales(tc.cort)[sindices]
  fbreaks,findices = freq_ticks(tc.cort.aspect,u[:,:,1])

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
                 time = vcat((ustrip.(eachindex(xi) * Δt(tempc))
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
