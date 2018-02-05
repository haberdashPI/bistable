using Sounds
using RCall
using Match
using Parameters
using ProgressMeter

struct CoherenceModel
  cort::CorticalModel
  ncomponents::Int
  window::Seconds{Float64}
  method::Any
  frame_len::Int
end

Δt(cohere::CoherenceModel) = Δt(cohere.cort)*cohere.frame_len
times(cohere::CoherenceModel,x) =
  times(cohere.cort,x)[min_window_size:cohere.frame_len:end]
times(cohere::CoherenceModel,C::EigenSeries) =
  (min_window_size:cohere.frame_len:length(C)*cohere.frame_len).*Δt(cohere.cort)
scales(cohere::CoherenceModel) = scales(cohere.cort)
rates(cohere::CoherenceModel) = rates(cohere.cort)

(cohere::CoherenceModel)(x::AbstractVector) = cohere(cohere.cort(x))
(cohere::CoherenceModel)(x::AbstractMatrix) = cohere(cohere.cort(x))

CoherenceModel(cort,ncomponents;window=1s,method=:pca,frame_len=10ms) =
  CoherenceModel(cort,ncomponents,window,method,
             max(1,floor(Int,frame_len/Δt(cort))))

const min_window_size = 10
windowing(x,dim;length=nothing,step=nothing) =
  (max(1,t-length):t for t in indices(x,dim)[min_window_size:step:end])

windowlen(cohere::CoherenceModel) = round(Int,cohere.window/Δt(cohere.cort))
nunits(cohere::CoherenceModel,x) = prod(size(x,3,4))
ncomponents(cohere::CoherenceModel) = cohere.ncomponents

# alternative: I could have a different
# set of eigenseries for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (cohere::CoherenceModel)(x)
  @match cohere.method begin
    :pca => begin
      windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                    step=cohere.frame_len))
      C = EigenSeries(eltype(x),length(windows),nunits(cohere,x),
                      ncomponents(cohere),Δt(cohere))

      @showprogress "Temporal Coherence Analysis: " for (i,w_inds) in windows
        window = x[w_inds,:,:,:]
        x_t = reshape(window,prod(size(window,1,2)),:)

        n = min(size(x_t,1),ncomponents(cohere))
        sv, = svds(x_t,nsv=n)

        λ = zeros(eltype(x),ncomponents(cohere))
        u = zeros(eltype(x),size(x_t,2),ncomponents(cohere))

        λ[1:n] = sv[:S].^2 ./ size(x_t,1)

        u[:,1:n] = sv[:V]
        var = mean(abs2.(x_t),1)
        C[i] = EigenSpace(sv[:V],(sv[:S]).^2 / size(x_t,1),var)
      end

      C
    end
    (:real_pca,n_phases) => begin
      windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                    step=cohere.frame_len))
      C = EigenSeries(real(eltype(x)),length(windows),nunits(cohere,x),
                      ncomponents(cohere),Δt(cohere))

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

        n = min(size(x_t,1),ncomponents(cohere))
        sv, = svds(x_t,nsv=n)

        λ = zeros(eltype(x),ncomponents(cohere))
        u = zeros(eltype(x),size(x_t,2),ncomponents(cohere))

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

# fusion_signal(cohere::CoherenceModel,C::EigenSeries) = fusion_signal(cohere.cort(x),C)
# function fusion_signal(cohere::CoherenceModel,C::EigenSeries)
#   vec(first.(eigvals.(C)) ./ sum.(var.(C)))
# end

const phase_resolution = 128

mask(cohere::CoherenceModel,C,x;kwds...) = mask(cohere,C,cohere.cort(x);kwds...)
"""
    mask(cohere,C,x;[phase=max_energy],[component=1])

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

* cohere - The settings for temporal coherence analysis.
* C - the previosuly computed temporal coherence, as an EigenSeries,
      (e.g. C = cohere(scene))
* x - the auditory scene
* phase - the phase of the component, see above
* component - the eigenvector of C to use (note that by default only the first
    is computed).
"""
function mask(cohere::CoherenceModel,C::EigenSpace,x::Array{T,4};
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

function __map_mask(fn::Function,cohere::CoherenceModel,
                    C::EigenSeries,x::Array{T},phase) where T
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                step=cohere.frame_len))
  y = zeros(length(windows))
  phases = zeros(y)

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m,phases[i] = mask(cohere,C[i],window,phase=phase)
    masked = inv(cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y,phases
end

function mask2(cohere::CoherenceModel,C::EigenSpace,x::Array{T,4};
               component=1) where T
  m = eigvecs(C)[:,component]
  m ./= maximum(abs.(m))
  m = reshape(m,size(x,3,4)...)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = m.*x[ii,:,:]
  end
  y
end

function __map_mask2(fn::Function,cohere::CoherenceModel,
                     C::EigenSeries,x::Array{T}) where T
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                step=cohere.frame_len))
  y = zeros(length(windows))

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask2(cohere,C[i],window)
    masked = inv(cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y
end

function mask(cohere::CoherenceModel,C::EigenSpace{<:Real},x::Array{T,4};
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

function __map_mask(fn::Function,cohere::CoherenceModel,
                    C::EigenSeries{<:Real},x::Array{T},
                    component) where T
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                step=cohere.frame_len))
  y = zeros(length(windows))

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask(cohere,C[i],window,component=component)
    masked = inv(cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y
end

function mean_spect(cohere::CoherenceModel,C::EigenSeries{<:Real},x::Array{T,4};
                    component=1) where T
  y = fill(real(zero(x[1])),size(x,1,4))
  norm = fill(zero(real(x[1])),size(x,1,4))
  dummy = __map_mask(cohere,C,x,component) do i,w_inds,masked,window
    y[w_inds,:] .+= masked
    norm[w_inds,:] .+= 1.0

    0.0
  end

  y ./ max.(1e-10,norm)
end

function scene_object_ratio(cohere::CoherenceModel,C::EigenSeries{<:Real},
                            x::Array{T,4},sp::Matrix;component=1) where T
  __map_mask(cohere,C,x,component) do i,w_inds,masked,window
    normed = sp[w_inds,:] ./ maximum(sp[w_inds,:])
    mean(masked .* normed) / mean(normed.^2)
  end
end

# TODO: make this work with mask2 (if we end up using that approach)
# function fusion_ratio2(cohere::CoherenceModel,C::EigenSeries,
#                       x::Array{T,4};phase=max_energy) where T

#   y,phases = __map_mask(cohere,C,x,phase) do i,w_inds,masked,window
#     sqrt(mean(abs2.(masked))) / sqrt(mean(abs2.(window)))
#   end
# end

function mean_spect(cohere::CoherenceModel,C::EigenSeries,x::Array{T,4}) where T
  y = fill(real(zero(x[1])),size(x,1,4))
  norm = fill(zero(real(x[1])),size(x,1,4))
  dummy = __map_mask2(cohere,C,x) do i,w_inds,masked,window
    y[w_inds,:] .+= masked
    norm[w_inds,:] .+= 1.0

    0.0
  end

  y ./ max.(1e-10,norm)
end
