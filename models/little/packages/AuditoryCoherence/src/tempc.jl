using Sounds
using RCall
using Match
using Parameters
using ProgressMeter
import AuditoryModel: raster_plot

struct CoherenceModel
  cort::CorticalModel
  ncomponents::Int
  window::Seconds{Float64}
  method::Any
  frame_len::Int
end

Δt(tc::CoherenceModel) = Δt(tc.cort)*tc.frame_len
times(tc::CoherenceModel,x) = times(tc.cort,x)[min_window_size:tc.frame_len:end]

(tc::CoherenceModel)(x::AbstractVector) = tc(tc.cort(x))
(tc::CoherenceModel)(x::AbstractMatrix) = tc(tc.cort(x))

CoherenceModel(cort,ncomponents;window=1s,method=:pca,frame_len=10ms) =
  CoherenceModel(cort,ncomponents,window,method,
             max(1,floor(Int,frame_len/Δt(cort))))

const min_window_size = 10
windowing(x,dim;length=nothing,step=nothing) =
  (max(1,t-length):t for t in indices(x,dim)[min_window_size:step:end])

windowlen(tc::CoherenceModel) = round(Int,tc.window/Δt(tc.cort))
nunits(tc::CoherenceModel,x) = prod(size(x,3,4))
ncomponents(tc::CoherenceModel) = tc.ncomponents

# alternative: I could have a different
# set of eigenseries for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (tc::CoherenceModel)(x)
  @match tc.method begin
    :pca => begin
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
    end
    (:real_pca,n_phases) => begin
      windows = enumerate(windowing(x,1;length=windowlen(tc),step=tc.frame_len))
      C = EigenSeries(real(eltype(x)),length(windows),nunits(tc,x),
                      ncomponents(tc),Δt(tc))

      @showprogress "Temporal Coherence Analysis: " for (i,w_inds) in windows
        window = x[w_inds,:,:,:]
        xc_t = reshape(window,prod(size(window,1,2)),:)

        # expand the representation, to have n_phase
        # real filters
        xp_t = Array{real(eltype(xc_t))}((size(xc_t,1),n_phases,size(xc_t,2)))
        for (j,phase) in enumerate(linspace(-π,π,n_phases+1)[1:end-1])
          xp_t[:,j,:] = real.(xc_t .* exp.(phase.*im))
        end
        x_t = reshape(xp_t,:,size(xc_t,2))

        n = min(size(x_t,1),ncomponents(tc))
        sv, = svds(x_t,nsv=n)

        λ = zeros(eltype(x),ncomponents(tc))
        u = zeros(eltype(x),size(x_t,2),ncomponents(tc))

        λ[1:n] = sv[:S].^2 ./ size(x_t,1)

        u[:,1:n] = sv[:V]
        var = squeeze(mean(abs2.(x_t),1),1)
        C[i] = EigenSpace(sv[:V],(sv[:S]).^2 / size(x_t,1),var)
      end

      C
    end
    method => begin
      error("No method named $(method)")
    end
  end
end

# fusion_signal(tc::CoherenceModel,C::EigenSeries) = fusion_signal(tc.cort(x),C)
# function fusion_signal(tc::CoherenceModel,C::EigenSeries)
#   vec(first.(eigvals.(C)) ./ sum.(var.(C)))
# end

const phase_resolution = 128

mask(tc::CoherenceModel,C,x;kwds...) = mask(tc,C,tc.cort(x);kwds...)
"""
    mask(tc,C,x;[phase=max_energy],[component=1])

Filter the given scene by a given principle component to extract
the object represented by that component from the scene.

Principle components are complex-valued, and so a given phase of the component
must be selected to find a real-valued mask. By default the maximum energy phase
is selected (by max_energy). You can also select the minimum energy phase
(min_energy), specifiy a specific phase value (e.g. π), or use a custom selector
function. The selector takes two arguments, `mask` and `phase` and should return
a utility function (higher is better), and will be used to find the maximum
utility phase. The mask is a vector representing the scale and frequency
responses in the same ordering as `vec(cr[1,1,:,:])` where `cr`
is a set of cortical responses.

# Arguments

* tc - The settings for temporal coherence analysis.
* C - the previosuly computed temporal coherence, as an EigenSeries,
      (e.g. C = tc(scene))
* x - the auditory scene
* phase - the phase of the component, see above
* component - the eigenvector of C to use (note that by default only the first
    is computed).
"""
function mask(tc::CoherenceModel,C::EigenSpace,x::Array{T,4};
              phase=max_energy,component=1) where T
  m,selected_phase = select_mask(C,x,phase,component)
  m = reshape(m,size(x,3,4)...)
  m ./= maximum(m)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = m.*x[ii,:,:]
  end
  y,selected_phase
end

select_mask(C::EigenSpace,x,phase::Number,component) =
  max.(0,real.(eigvecs(C)[:,component] .* exp.(phase*im))), phase

function select_mask(C::EigenSpace{T},x,phase_selector,component) where T
  phases = linspace(-π,π,phase_resolution+1)[1:end-1]
  vals = similar(phases)
  @showprogress "evaluating phases: " for (i,p) in enumerate(phases)
    vals[i] = phase_selector(x,select_mask(C,x,p,component)...)
  end
  i = indmax(vals)
  select_mask(C,x,phases[i],component)
end

max_energy(x,mask,phase) = sum(mask.^2)
min_energy(x,mask,phase) = -sum(mask.^2)

function max_filtering(x,mask,phase)
  mask = reshape(mask,size(x,3,4)...)
  mask = mask ./ maximum(mask)
  sum(CartesianRange(size(x,1,2))) do ii
    sum(abs2.(mask.*x[ii,:,:]))
  end
end
min_filtering(x,mask,phase) = -max_filtering(x,mask,phase)

function __map_mask(fn::Function,tc::CoherenceModel,
                    C::EigenSeries,x::Array{T},phase) where T
  windows = enumerate(windowing(x,1;length=windowlen(tc),step=tc.frame_len))
  y = zeros(length(windows))
  phases = zeros(y)

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m,phases[i] = mask(tc,C[i],window,phase=phase)
    masked = inv(cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y,phases
end

function mask2(tc::CoherenceModel,C::EigenSpace,x::Array{T,4};component=1) where T
  m = eigvecs(C)[:,component]
  m ./= maximum(abs.(m))
  m = reshape(m,size(x,3,4)...)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = m.*x[ii,:,:]
  end
  y
end

function __map_mask2(fn::Function,tc::CoherenceModel,
                     C::EigenSeries,x::Array{T}) where T
  windows = enumerate(windowing(x,1;length=windowlen(tc),step=tc.frame_len))
  y = zeros(length(windows))

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask2(tc,C[i],window)
    masked = inv(cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y
end

function mask(tc::CoherenceModel,C::EigenSpace{<:Real},x::Array{T,4};
              component=1) where T
  m = if component == :max
    eigvecs(C)[:,indmax(eigvals(C))]
  else
    eigvecs(C)[:,component]
  end
  m ./= maximum(abs.(m))
  m .= max.(0,m)
  mr = reshape(m,size(x,3,4)...)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = mr.*x[ii,:,:]
  end
  y
end

function __map_mask(fn::Function,tc::CoherenceModel,
                    C::EigenSeries{<:Real},x::Array{T},
                    component) where T
  windows = enumerate(windowing(x,1;length=windowlen(tc),step=tc.frame_len))
  y = zeros(length(windows))

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask(tc,C[i],window,component=component)
    masked = inv(cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y
end

function mean_spect(tc::CoherenceModel,C::EigenSeries{<:Real},x::Array{T,4};
                    component=1) where T
  y = fill(real(zero(x[1])),size(x,1,4))
  norm = fill(zero(real(x[1])),size(x,1,4))
  dummy = __map_mask(tc,C,x,component) do i,w_inds,masked,window
    y[w_inds,:] .+= masked
    norm[w_inds,:] .+= 1.0

    0.0
  end

  y ./ max.(1e-10,norm)
end

function scene_object_ratio(tc::CoherenceModel,C::EigenSeries{<:Real},
                            x::Array{T,4},sp::Matrix;component=1) where T
  __map_mask(tc,C,x,component) do i,w_inds,masked,window
    normed = sp[w_inds,:] ./ maximum(sp[w_inds,:])
    mean(masked .* normed) / mean(normed.^2)
  end
end

# TODO: make this work with mask2 (if we end up using that approach)
# function fusion_ratio2(tc::CoherenceModel,C::EigenSeries,
#                       x::Array{T,4};phase=max_energy) where T

#   y,phases = __map_mask(tc,C,x,phase) do i,w_inds,masked,window
#     sqrt(mean(abs2.(masked))) / sqrt(mean(abs2.(window)))
#   end
# end

function mean_spect(tc::CoherenceModel,C::EigenSeries,x::Array{T,4}) where T
  y = fill(real(zero(x[1])),size(x,1,4))
  norm = fill(zero(real(x[1])),size(x,1,4))
  dummy = __map_mask2(tc,C,x) do i,w_inds,masked,window
    y[w_inds,:] .+= masked
    norm[w_inds,:] .+= 1.0

    0.0
  end

  y ./ max.(1e-10,norm)
end

rplot(tc::CoherenceModel,x::Sound;kwds...) = rplot(tc,tc(x);kwds...)

function rplot(tc::CoherenceModel,λ::Vector)
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
function rplot(tc::CoherenceModel,C::EigenSeries;n=ncomponents(C))
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

function rplot(tc::CoherenceModel,C::EigenSpace;
               n=ncomponents(C),showvar=true,λ_digits=:automatic)
  λ = abs.(eigvals(C))
  order = sortperm(λ,rev=true)
  λ = λ[order]
  u = eigvecs(C)[:,order]
  u = u[:,1:min(n,end)]
  if showvar; u = [u var(C)]; end
  u = reshape(u,length(scales(tc.cort)),:,size(u,2))
  ii = CartesianRange(size(u))
  at(i) = vec(map(ii -> ii[i],ii))
  digits = λ_digits == :automatic ? -floor(Int,log10(λ[end]))+2 : λ_digits

  function title(n)
    if n <= length(λ)
      # nstr = @sprintf("%02d",n)
      "lambda[$(n)] == $(round(λ[n],digits))"
      # "Lmb_$nstr = $(round(λ[n],3))"
    else
      "sigma^2 == $(round(sum(var(C)),digits))"
      # "Variance = $(sum(var(C)))"
    end
  end

  df = DataFrame(response = vec(u),
                 scale_index = at(1),
                 freq_bin = at(2),
                 component = at(3),
                 component_title = title.(at(3)))

  sindices = 1:2:length(scales(tc.cort))
  sbreaks = round.(scales(tc.cort)[sindices],2)
  fbreaks,findices = freq_ticks(tc.cort.aspect,u[:,:,1])

  p = raster_plot(df,value=:response,x=:scale_index,y=:freq_bin)

R"""

  library(ggplot2)

  $p + facet_wrap(~component_title,labeller=label_parsed) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    scale_x_continuous(breaks=$sindices,labels=$sbreaks) +
    ylab('Frequency (Hz)') + xlab('Scale (cycles/octave)')
"""

end

function rplot(tempc::CoherenceModel,C::Array{<:EigenSeries},
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
