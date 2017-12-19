# TODO: for updating via adaptation I think it should be possible to think in
# terms of the updates happening to the PCA representation, and then show the
# equivalence of that to the updates happening to the covariance matrix
# (and/or, do we care if it is equivalent, perhaps this is the representation
# we want to operate over)

using Parameters
using ProgressMeter
include("cortical.jl")
include("online_pca.jl")

abstract type TCAnalysis end
Δt(tc::TCAnalysis) = Δt(tc.upstream)
times(tc::TCAnalysis,x) = times(tc.upstream,x)

struct OnlineTCAnalysis <: TCAnalysis
  upstream::CorticalModel
  ncomponents::Int
  rate::Seconds{Float64}
  method::Symbol
end

(tc::OnlineTCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
(tc::OnlineTCAnalysis)(x::AbstractMatrix) = tc(tc.upstream(x))

TCAnalysis(upstream,ncomponents,rate=1s;method=:pca) =
  OnlineTCAnalysis(upstream,ncomponents,rate,method)

# alternative: I could have a different
# set of eigenseries for each time scale
Base.CartesianRange(x::Int) = CartesianRange((x,))
function (tc::OnlineTCAnalysis)(x)
  if tc.method == :pca
    C = EigenSeries(eltype(x),size(x,1),prod(size(x,3,4)),tc.ncomponents)
    window_len = round(Int,tc.rate/Δt(tc))
    @showprogress for t in indices(x,1)
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
      var = sum(abs2.(x_t)) / size(x_t,1)
      C[t] = EigenSpace(sv[:V],(sv[:S]).^2 / size(x_t,1),var)
    end

    C
  elseif tc.method == :ipca
    C_t = EigenSpace(eltype(x),prod(size(x,3,4)),tc.ncomponents)
    C = EigenSeries(size(x,1),C_t)

    dt = Δt(tc) / tc.rate

    function helper__(C_t,C,x,dt)

      @showprogress for t in indices(x,1)
        for r in indices(x,2)
          # approximately: C_t = (1-dt)C_t + x*x'*dt
          C_t = update(C_t,x[t,r,:,:],dt)
        end
        C[t] = C_t
      end

      C
    end

    helper__(C_t,C,x,dt)
  elseif tc.method == :ipca_lr
    C_t = EigenSpace(eltype(x),prod(size(x,3,4)),tc.ncomponents)
    C = EigenSeries(size(x,1),C_t)

    dt = Δt(tc) / tc.rate

    function helper___(C_t,C,x,dt)

      @showprogress for t in indices(x,1)
        for r in indices(x,2)
          # approximately: C_t = (1-dt)C_t + x*x'*dt
          C_t = update_approx(C_t,x[t,r,:,:],dt)
        end
        C[t] = C_t
      end

      C
    end

    helper___(C_t,C,x,dt)
  elseif tc.method == :ipca_rsplit
    nr = length(rates(tc.upstream))
    C_t = [EigenSpace(eltype(x),prod(size(x,3,4)),tc.ncomponents) for r in 1:nr]
    C = [EigenSeries(size(x,1),C_t[r]) for r in 1:nr]
    dt = Δt(tc) / tc.rate

    function helper_(C_t,C,x,dt)
      @showprogress for t in indices(x,1)
        for r in indices(x,2)
          # approximately: C_t = (1-dt)C_t + x*x'*dt
          C_t[r] = update(C_t[r],x[t,r,:,:],dt)
          C[r][t] = C_t[r]
        end
      end
      C
    end

    helper_(C_t,C,x,dt)
  else
    error("No method named $(tc.method)")
  end
end

fusion_signal(tc::OnlineTCAnalysis,C::EigenSeries,x::AbstractVector) =
  fusion_signal(tc.upstream(x),C,x)
function fusion_signal(tc::OnlineTCAnalysis,C::EigenSeries,x)
  vec(first.(eigvals.(C)) ./ var.(C))
  # cumvar = cumsum(sum(abs2.(x),[2,3,4]),1) ./ (size(x,2)*indices(x,1))
  # vec(first.(eigvals.(C)) ./ cumvar)
end

rplot(tc::OnlineTCAnalysis,x::TimedSound.Sound;kwds...) = rplot(tc,tc(x);kwds...)

function rplot(tc::OnlineTCAnalysis,λ::Vector)
  @assert all(imag.(λ) .== 0) "Can't plot complex eigenvalues"
  λ = real.(λ)
  df = DataFrame(value = sort(λ,rev=true),index = collect(eachindex(λ)))

R"""
  library(ggplot2)

  ggplot($df,aes(x=index,y=value)) + geom_bar(stat='identity')
"""
end

function rplot(tc::OnlineTCAnalysis,C::EigenSeries;n=ncomponents(C),
               oddonly=false)
  λ = C.λ
  λ = λ[:,sortperm(abs.(λ[max(1,end-10),:]),rev=true)]
  ii = CartesianRange(size(λ))
  at(i) = map(ii -> ii[i],ii)

  prop_complex = mean(imag.(λ))
  if prop_complex > 0
    warn("$(round(100prop_complex,1))% of eigenvalues were "*
         "complex (using absolute value).")
  end

  df = DataFrame(value = vec(abs.(λ)),
                 time = vec(ustrip(at(1) * Δt(spect))),
                 component = vec(at(2)))

  if oddonly
    df = df[isodd.(df[:component]),:]
  end
  df = df[df[:component] .<= n,:]

R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=value,color=factor(component),group=component)) +
    geom_line() + scale_color_brewer(palette='Set1',name='Component') +
    xlab('Time (s)') + ylab('Value')
"""
end



function rplot(tc::OnlineTCAnalysis,C::EigenSpace;n=ncomponents(C),oddonly=false)
  λ = abs.(eigvals(C))
  order = sortperm(λ,rev=true)
  λ = λ[order]
  u = eigvecs(C)[:,order]
  u = u[:,1:min(n,end)]
  u = reshape(u,length(scales(tc.upstream)),:,size(u,2))
  ii = CartesianRange(size(u))
  at(i) = vec(map(ii -> ii[i],ii))
  function title(n)
    nstr = @sprintf("%02d",n)
    "Lmb_$nstr = $(round(λ[n],3))"
  end

  colormap = "#".*hex.(RGB.(cmap("C1")))

  df = DataFrame(r_phase = angle.(vec(u)),
                 r_amp = abs.(vec(u)),
                 scale_index = at(1),
                 freq_bin = at(2),
                 component = title.(at(3)))

  fbreaks = 2.0.^(-3:2)
  fs = freqs(tc.upstream,u[:,:,1])
  findices = mapslices(abs.(1000.0.*fbreaks .- fs'),2) do row
    _, i = findmin(row)
    i
  end

  if oddonly
    df = df[isodd.(df[:component]),:]
  end

  sindices = 1:2:length(scales(tc.upstream))
  sbreaks = scales(tc.upstream)[sindices]

R"""

  library(ggplot2)

  ggplot($df,aes(x=scale_index,y=freq_bin,fill=r_phase,alpha=r_amp)) +
    geom_raster() + facet_wrap(~component) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    scale_x_continuous(breaks=$sindices,labels=$sbreaks) +
    ylab('Frequency (kHz)') + xlab('Scale') +
    scale_fill_gradientn(colors=$colormap) +
    scale_alpha_continuous(range=c(0,1))

"""

end

# struct BatchTCAnalysis <: TCAnalysis
#   upstream::CorticalModel
#   ncomponents::Int
# end

# (tc::BatchTCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
# function (tc::BatchTCAnalysis)(x)
#   sv, = svds(reshape(x,prod(size(x,1,2)),:),nsv=tc.ncomponents)

#   yv = reshape(x,prod(size(x,1,2)),:)*sv[:V]
#   y = reshape(yv,size(x,1:2...)...,:)
#   return sv[:S].^2 / size(x,1), reshape(sv[:V]',:,size(x,3,4)...), y
# end

# TCAnalysis(upstream,ncomponents;method=:ipca,init_len=0) =
#   method ∈ [:batch,:pca] ? BatchTCAnalysis(upstream,ncomponents) :
#   OnlineTCAnalysis(upstream,ncomponents,method,init_len)

#========================================
struct TCAnalysis
  upstream::CorticalModel
  ncomponents::Int
end
times(tc::TCAnalysis,x) = times(tc.upstream,x)
times(tc::TCAnalysis,x::Array{T,3}) where T = times(tc.upstream,x[:,1,:])
freqs(tc::TCAnalysis,x::Array{T,3}) where T = freqs(tc.upstream,x[:,1,:])

(tc::TCAnalysis)(x::AbstractVector) = tc(tc.upstream(x))
function (tc::TCAnalysis)(x)
  sv, = svds(reshape(x,prod(size(x,1,2)),:),nsv=tc.ncomponents)

  sv[:S].^2 / size(x,1), reshape(sv[:V],:,size(x,3,4))
end
========================================#

# function rplot(tc::TCAnalysis,data::Matrix)
#   ixs = CartesianRange(size(abs.(data)))
#   at(ixs,i) = map(x -> x[i],ixs)

#   df = DataFrame(response = vec(data),
#                  time = vec(times(tc,data)[at(ixs,1)]),
#                  component_index = vec(at(ixs,2)))

# R"""

#   library(ggplot2)

#   ggplot($df,aes(x=time,y=component_index,fill=response)) +
#     geom_raster() +
#     ylab('Component') + ('Time (s)') +
#     scale_fill_distiller(palette='Reds',name='Amplitude',
#                          direction=1)

# """
# end

# function rplot(tc::TCAnalysis,data::Array{T,3}) where T
#   ixs = CartesianRange(size(abs.(data)))
#   at(ixs,i) = map(x -> x[i],ixs)

#   df = DataFrame(response = vec(data),
#                  time = vec(times(tc,data)[at(ixs,1)]),
#                  component_index = vec(at(ixs,2)),
#                  freq_index = vec(at(ixs,3)))

#   fbreaks = 2.0.^(-3:2)
#   fs = freqs(tc,data)
#   findices = mapslices(abs.(1000.0.*fbreaks .- fs'),2) do row
#     _, i = findmin(row)
#     i
#   end

# R"""

#   library(ggplot2)

#   ggplot($df,aes(x=time,y=component_index,fill=response)) +
#     geom_raster() +
#     facet_wrap(~component_index) +
#     scale_y_continuous(breaks=$findices,labels=$fbreaks) +
#     ylab('Frequency (kHz') + xlab('Time (s)') +
#     scale_fill_distiller(palette='Reds',name='Amplitude',
#                          direction=1)

# """
# end
