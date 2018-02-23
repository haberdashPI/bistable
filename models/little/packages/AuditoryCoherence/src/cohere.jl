using Sounds
using RCall
using Match
using Parameters
using ProgressMeter
using Parameters

export NMFC, NMFDirect, realpos

realpos(x) = max(0,real(x))

abstract type CoherenceModelMethod end
struct CoherenceModel{T} <: CoherenceModelMethod
  cort::CorticalModel
  ncomponents::Int
  window::Seconds{Float64}
  minwindow::Seconds{Float64}
  method::T
  delta::Seconds{Float64}
end

@with_kw struct CoherencePCA
  normalize::Bool = true
end
@with_kw struct CoherenceRealPCA
  n_phases::Int = 12
  normalize::Bool = true
end

struct NMFDirect end
@with_kw struct NMFC
  nonlinear = realpos
end

@with_kw struct CoherenceNMF{M}
  submethod::M = NMFDirect()
  normalize::Bool = true
  maxiter::Int = 2000
  tol::Float64 = 1e-4
  normalize_tc::Seconds{Float64} = 1s
end
nmf_tc(cohere::CoherenceModel{CoherenceNMF{T}}) where T =
  1 / floor(Int,cohere.method.normalize_tc/Δt(cohere.cort))

Δt(cohere::CoherenceModel) = cohere.delta
frame_length(cohere::CoherenceModel) =
  min(1,floor(Int,cohere.delta / Δt(cohere.cort)))
times(cohere::CoherenceModel,x::AbstractArray) =
  times(cohere.cort,x)[min_windowlen(cohere):frame_length(cohere):end]
times(cohere::CoherenceModel,C::FactorSeries) =
  ((0:length(C)-1).*frame_length(cohere) .+ min_windowlen(cohere)) .*
  Δt(cohere.cort)
scales(cohere::CoherenceModel) = scales(cohere.cort)
rates(cohere::CoherenceModel) = rates(cohere.cort)

(cohere::CoherenceModel)(x::AbstractVector) = cohere(cohere.cort(x))
(cohere::CoherenceModel)(x::AbstractMatrix) = cohere(cohere.cort(x))

function CoherenceModel(cort,ncomponents;window=1s,minwindow=window,
                        method=:pca,delta=10ms,
                        normalize_phase=true,method_kwds...)

  method = @match method begin
    :pca => CoherencePCA(;method_kwds...)
    :real_pca => CoherenceRealPCA(;method_kwds...)
    :nmf => CoherenceNMF(;method_kwds...)
  end

  CoherenceModel(cort,ncomponents,convert(Seconds{Float64},window),
                 convert(Seconds{Float64},minwindow),method,
                 convert(Seconds{Float64},delta))
end

windowing(x,dim;length=nothing,step=nothing,minlength=nothing) =
  (max(1,t-length+1):t for t in indices(x,dim)[minlength:step:end])

windowlen(cohere::CoherenceModel) = round(Int,cohere.window/Δt(cohere.cort))
min_windowlen(cohere::CoherenceModel) =
  round(Int,cohere.minwindow / Δt(cohere.cort))
nunits(cohere::CoherenceModel,x) = prod(size(x,3,4))
ncomponents(cohere::CoherenceModel) = cohere.ncomponents

# alternative: I could have a different
# set of eigenseries for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (cohere::CoherenceModel{CoherencePCA})(x)
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                minlength=min_windowlen(cohere),
                                step=frame_length(cohere)))
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

  mehtod.normalize ? normalize_phase!(C) : C
end

function (cohere::CoherenceModel{CoherenceRealPCA})(x)
  n_phases = cohere.method.n_phases
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                minlength=min_windowlen(cohere),
                                step=frame_length(cohere)))
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

  method.normalize ? normalize_phase!(C) : C
end

function (cohere::CoherenceModel{CoherenceNMF{NMFDirect}})(x)
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                minlength=min_windowlen(cohere),
                                step=frame_length(cohere)))

  C = NMFSeries(length(windows),windowlen(cohere)*length(rates(cohere)),
                nunits(cohere,x),ncomponents(cohere),Δt(cohere.cort),Δt(cohere))
  convergence_count = 0
  # Winit,Hinit = fill(0.0,(0,0)),fill(0.0,(0,0))
  # method = NMF.ALSPGrad{Float64}(tol=cohere.method.tol,
  #                                maxiter=cohere.method.maxiter)

  @showprogress "Temporal Coherence Analysis: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    x_t = reshape(window,prod(size(window,1,2)),:)

    k = min(size(x_t,1),ncomponents(cohere))

    # if isempty(Winit)
    #   Winit,Hinit = NMF.nndsvd(abs.(x_t),k)
    # end

    # solution = NMF.solve!(method,abs.(x_t),Winit,Hinit)
    solution = nnmf(abs.(x_t) .+ cohere.method.tol/2,k,init=:nndsvd,
                    tol=cohere.method.tol,
                    maxiter=cohere.method.maxiter)

    convergence_count += !solution.converged
    solution.converged || warn("NMF failed to converge.",once=true,
                               key=object_id(C))

    # Winit,Hinit = solution.W,solution.H
    C[i] = NMFSpace(solution.W,solution.H,Δt(cohere.cort))
  end

  if convergence_count > 0
    info("$(100round(convergence_count / length(windows),3))% of frames "*
         " failed to fully converge to a solution.")
  end

  cohere.method.normalize ? normalize_components!(C,nmf_tc(cohere)) : C
end

function (cohere::CoherenceModel{CoherenceNMF{NMFC}})(x)
  C = NMFSeries(size(x,1),nunits(cohere,x)*length(rates(cohere)),
                nunits(cohere,x),ncomponents(cohere),Δt(cohere.cort),
                Δt(cohere.cort))
  convergence_count = 0
  Winit,Hinit = fill(0.0,(0,0)),fill(0.0,(0,0))
  nonlinear = cohere.method.submethod.nonlinear

  @showprogress "Temporal Coherence Analysis: " for t in 1:size(x,1)
    x_t = reshape(x[t,:,:,:],size(x,2),:)
    Ci = Array{eltype(x_t)}(size(x,2),nunits(cohere,x),nunits(cohere,x))
    for r in 1:size(x_t,1) Ci[r,:,:] = x_t[r,:] .* x_t[r,:]' end
    Cir = reshape(Ci,:,nunits(cohere,x))

    k = min(size(Ci,1),ncomponents(cohere))

    if isempty(Winit)
      Winit,Hinit = NMF.nndsvd(nonlinear.(Cir),k)
    end

    solution = nnmf(nonlinear.(Cir),k;tol=cohere.method.tol,
                    maxiter=cohere.method.maxiter)

    convergence_count += !solution.converged
    solution.converged || warn("NMF failed to converge.",once=true,
                               key=object_id(C))

    Winit,Hinit = solution.W,solution.H
    C[t] = NMFSpace(solution.W,solution.H,Δt(cohere.cort))
  end

  if convergence_count > 0
    info("$(100round(convergence_count / size(x,1),3))% of frames "*
         " failed to fully converge to a solution.")
  end

  cohere.method.normalize ? normalize_components(C,nmf_tc(cohere)) : C
end

function mask(cohere::CoherenceModel,C::NMFSpace,
              x::AbstractArray{T,4} where T;component=1,maxvalue=maximum(abs,C))
  @assert nunits(C) == prod(size(x,3,4))
  y = copy(x)
  c_ = factors(C)[:,component]
  c = reshape(c_,length(scales(cohere)),:)
  c ./= maxvalue
  @simd for ii in CartesianRange(size(x,1,2))
    @inbounds y[ii,:,:] .= sqrt.(abs.(y[ii,:,:]) .* c) .*
      exp.(angle.(y[ii,:,:])*im)
  end
  y
end

function mask(cohere::CoherenceModel,C::EigenSpace{<:Complex},
              x::AbstractArray{T,4} where T;component=1,maxvalue=maximum(abs,C))
  @assert nunits(C) == prod(size(x,3,4))
  y = copy(x)
  pc_ = eigvecs(C)[:,component]
  pc = reshape(pc_,length(scales(cohere)),:)
  pc ./= maxvalue
  @simd for ii in CartesianRange(size(x,1,2))
    @inbounds y[ii,:,:] .= sqrt.(abs.(y[ii,:,:]) .* max.(0,real.(pc))) .*
      exp.(angle.(y[ii,:,:])*im)
  end
  y
end

function mask(cohere::CoherenceModel,C::EigenSpace{<:Real},x::AbstractArray{T,4};
              component=1,maxvalue=maximum(abs,C)) where T
  m = if component == :max
    eigvecs(C)[:,indmax(eigvals(C))]
  else
    eigvecs(C)[:,component]
  end
  m ./= maxvalue
  m .= max.(0,m)
  mr = reshape(m,size(x,3,4)...)

  y = similar(x)
  for ii in CartesianRange(size(x,1,2))
    y[ii,:,:] = mr.*x[ii,:,:]
  end
  y
end

function __map_mask(fn::Function,cohere::CoherenceModel,C::FactorSeries,
                    x::AbstractArray{T} where T,component)
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                minlength=min_windowlen(cohere),
                                step=frame_length(cohere)))
  y = zeros(length(windows))
  max_C = maximum(abs,C)
  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask(cohere,C[i],window,component=component,maxvalue=max_C)
    masked = inv(cohere.cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y
end

function mean_spect2(cohere::CoherenceModel,C::NMFSeries,
                     x::AbstractArray{T,4};component=1) where T
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                minlength=min_windowlen(cohere),
                                step=frame_length(cohere)))
  y = fill(zero(x[1]),size(x))
  norm = fill(real(zero(x[1])),size(x))
  @showprogress "Masking: " for (i,w_inds) in windows
    c = factors(C[i])[:,component]
    y[w_inds,:,:,:] .+= reshape(c,1,1,length(scales(cohere)),:)
    norm[w_inds,:,:,:] += 1
  end
  y ./= norm
  y ./= maximum(abs,y)
  y .= sqrt.(abs.(x) .* y) .* exp.(angle.(x).*im)

  inv(cohere.cort,y)
end

function mean_spect(cohere::CoherenceModel,C::FactorSeries,
                    x::AbstractArray{T,4};component=1) where T
  y = fill(real(zero(x[1])),size(x,1,4))
  norm = fill(zero(real(x[1])),size(x,1,4))
  dummy = __map_mask(cohere,C,x,component) do i,w_inds,masked,window
    y[w_inds,:] .+= masked
    norm[w_inds,:] .+= 1.0

    0.0
  end

  y ./ max.(1e-10,norm)
end
