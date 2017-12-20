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
  elseif tc.method == :ipca
    C_t = EigenSpace(eltype(x),prod(size(x,3,4)),tc.ncomponents)
    C = EigenSeries(size(x,1),C_t,Δt(tc))

    dt = Δt(tc) / tc.rate

    function helper__(C_t,C,x,dt)

      @showprogress "Temporal Coherence Analysis: " for t in indices(x,1)
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
    C = EigenSeries(size(x,1),C_t,Δt(tc))

    dt = Δt(tc) / tc.rate

    function helper___(C_t,C,x,dt)

      @showprogress "Temporal Coherence Analysis: " for t in indices(x,1)
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
    C = [EigenSeries(size(x,1),C_t[r],Δt(tc)) for r in 1:nr]
    dt = Δt(tc) / tc.rate

    function helper_(C_t,C,x,dt)
      @showprogress "Temporal Coherence Analysis: " for t in indices(x,1)
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

fusion_signal(tc::TCAnalysis,C::EigenSeries,x::AbstractVector) =
  fusion_signal(tc.upstream(x),C,x)
function fusion_signal(tc::TCAnalysis,C::EigenSeries,x)
  vec(first.(eigvals.(C)) ./ sum.(var.(C)))
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

function rplot(tc::TCAnalysis,C::EigenSeries;n=ncomponents(C),
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

function rplot(tc::TCAnalysis,C::EigenSpace;n=ncomponents(C),oddonly=false)
  λ = abs.(eigvals(C))
  order = sortperm(λ,rev=true)
  λ = λ[order]
  u = eigvecs(C)[:,order]
  u = u[:,1:min(n,end)]
  u = [u var(C)]
  u = reshape(u,length(scales(tc.upstream)),:,size(u,2))
  ii = CartesianRange(size(u))
  at(i) = vec(map(ii -> ii[i],ii))
  function title(n)
    if n <= length(λ)
      nstr = @sprintf("%02d",n)
      "Lmb_$nstr = $(round(λ[n],3))"
    else
      "Variance = $(sum(var(C)))"
    end
  end

  colormap = "#".*hex.(RGB.(cmap("C1")))

  df = DataFrame(r_phase = angle.(vec(u)),
                 r_amp = abs.(vec(u)),
                 scale_index = at(1),
                 freq_bin = at(2),
                 component = title.(at(3)))

  if oddonly
    df = df[isodd.(df[:component]),:]
  end

  sindices = 1:2:length(scales(tc.upstream))
  sbreaks = scales(tc.upstream)[sindices]
  fbreaks,findices = freq_ticks(tc.upstream.aspect,u[:,:,1])

R"""

  library(ggplot2)

  ggplot($df,aes(x=scale_index,y=freq_bin,fill=r_phase,alpha=r_amp)) +
    geom_raster() + facet_wrap(~component) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    scale_x_continuous(breaks=$sindices,labels=$sbreaks) +
    ylab('Frequency (Hz)') + xlab('Scale') +
    scale_fill_gradientn(colors=$colormap) +
    scale_alpha_continuous(range=c(0,1))

"""

end

function plot_resps(x,vars,varname)
  df = DataFrame(resp = vcat((real.(x[i])
                              for i in 1:length(x))...),
                 var = vcat((fill(vars[i],length(x[i]))
                             for i in eachindex(x))...),
                 time = vcat((ustrip.(eachindex(xi) * Δt(spect))
                              for xi in x)...))

R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=resp,group=factor(var),color=factor(var))) +
    geom_line() +
    scale_color_brewer(palette='Set1',name=$varname) +
    coord_cartesian(ylim=c(1.1,0)) + ylab('lambda / var(x)') +
    xlab('time (s)')
"""
end
