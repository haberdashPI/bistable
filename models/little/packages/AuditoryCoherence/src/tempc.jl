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
  normalize::Bool
end

const min_window_size = 10

Δt(cohere::CoherenceModel) = Δt(cohere.cort)*cohere.frame_len
function times(cohere::CoherenceModel,x)
  @show length(times(cohere.cort,x))
  times(cohere.cort,x)[min_window_size:cohere.frame_len:end]
end
times(cohere::CoherenceModel,C::EigenSeries) =
  ((0:length(C)-1).*cohere.frame_len .+ min_window_size) .* Δt(cohere.cort)
scales(cohere::CoherenceModel) = scales(cohere.cort)
rates(cohere::CoherenceModel) = rates(cohere.cort)

(cohere::CoherenceModel)(x::AbstractVector) = cohere(cohere.cort(x))
(cohere::CoherenceModel)(x::AbstractMatrix) = cohere(cohere.cort(x))

CoherenceModel(cort,ncomponents;window=1s,method=:pca,frame_len=10ms,
               normalize_phase=true) =
  CoherenceModel(cort,ncomponents,window,method,
                 max(1,floor(Int,frame_len/Δt(cort))),normalize_phase)

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

      cohere.normalize ? normalize_phase!(C) : C
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

      cohere.normalize ? normalize_phase!(C) : C
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

function mask(cohere::CoherenceModel,C::EigenSpace{<:Complex},
              x::AbstractArray{T,4} where T;component=1)
  @assert nfeatures(C) == prod(size(x,3,4))
  y = copy(x)
  pc_ = eigvecs(C)[:,component]
  pc = reshape(pc_,length(scales(cohere)),:)
  pc ./= maximum(abs.(pc))
  @simd for ii in CartesianRange(size(x,1,2))
    @inbounds y[ii,:,:] .= sqrt.(abs.(y[ii,:,:]) .* max.(0,real.(pc))) .*
      exp.(angle.(y[ii,:,:])*im)
  end
  y
end

function mask(cohere::CoherenceModel,C::EigenSpace{<:Real},x::AbstractArray{T,4};
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

function __map_mask(fn::Function,cohere::CoherenceModel,C::EigenSeries,
                    x::AbstractArray{T} where T,component)
  windows = enumerate(windowing(x,1;length=windowlen(cohere),
                                step=cohere.frame_len))
  y = zeros(length(windows))

  @showprogress "Calculating Mask: " for (i,w_inds) in windows
    window = x[w_inds,:,:,:]
    m = mask(cohere,C[i],window,component=component)
    masked = inv(cohere.cort,m)

    masked ./= maximum(abs.(masked))
    window ./= maximum(abs.(window))
    y[i] = fn(i,w_inds,masked,window)
  end

  y
end

function mean_spect(cohere::CoherenceModel,C::EigenSeries,x::AbstractArray{T,4};
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
