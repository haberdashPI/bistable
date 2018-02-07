using RCall
import AuditoryModel: rplot, raster_plot
export rplot, scale_plot
R"library(ggplot2)"

rplot(cohere::CoherenceModel,x::Sound;kwds...) = rplot(cohere,cohere(x);kwds...)

function rplot(cohere::CoherenceModel,λ::Vector)
  @assert all(imag.(λ) .== 0) "Can't plot complex eigenvalues"
  λ = real.(λ)
  df = DataFrame(value = sort(λ,rev=true),index = collect(eachindex(λ)))

R"""
  ggplot($df,aes(x=index,y=value)) + geom_bar(stat='identity')
"""
end

# TODO: improve this plot by making it relative to variance (ala fusion signal)
# and then make plot_resps a version of this using an array of eigenseries
function rplot(cohere::CoherenceModel,C::EigenSeries;n=ncomponents(C))
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
                 time = vec(ustrip(at(1) * Δt(cohere))),
                 component = vec(at(2)))

  df = df[df[:component] .<= n,:]

R"""
  ggplot($df,aes(x=time,y=value,color=factor(component),group=component)) +
    geom_line() + scale_color_brewer(palette='Set1',name='Component') +
    xlab('Time (s)') + ylab('Value')
"""
end

function scale_plot(cohere::CoherenceModel,C::EigenSeries;
                    components=1:ncomponents(C),scales=AuditoryModel.scales(cohere))
  nscales = length(AuditoryModel.scales(cohere))
  scales,indices = AuditoryModel.findnear(AuditoryModel.scales(cohere),scales)

  u = reshape(eigvecs(C),length(C),nscales,:,ncomponents(C))
  u = u[:,indices,:,components]
  ii = CartesianRange(size(u))
  at(i) = map(ii -> ii[i],ii)

  df = DataFrame(value = vec(u),
                 time = vec(ustrip.(times(cohere,C)[at(1)])),
                 scale = vec(round.(scales[at(2)],2)),
                 freq_bin = vec(at(3)),
                 component = vec(at(4)))

  fbreaks,findices = freq_ticks(cohere.cort.aspect)
  p = raster_plot(df,value=:value,x=:time,y=:freq_bin)

R"""

  scalestr = function(x){sprintf("Scale: %3.2f cyc/oct",x)}
  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr($scales))
  }

  $p + facet_grid(ordered_scales(scale)~paste("component",component)) +
  scale_y_continuous(breaks=$findices,labels=$fbreaks) +
  ylab('Frequency (Hz)') + xlab('Time (s)')

"""
end

function rplot(cohere::CoherenceModel,C::EigenSpace;
               n=ncomponents(C),showvar=true,λ_digits=:automatic)
  λ = abs.(eigvals(C))
  order = sortperm(λ,rev=true)
  λ = λ[order]
  u = eigvecs(C)[:,order]
  u = u[:,1:min(n,end)]
  if showvar; u = [u var(C)]; end
  u = reshape(u,length(scales(cohere.cort)),:,size(u,2))
  ii = CartesianRange(size(u))
  at(i) = vec(map(ii -> ii[i],ii))
  digits = λ_digits != :automatic ?  λ_digits :
    min(-floor(Int,log10(λ[minimum(components)]))+2,
        -floor(Int,log10(λ[maximum(components)]))+1)

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

  sindices = 1:2:length(scales(cohere.cort))
  sbreaks = round.(scales(cohere.cort)[sindices],2)
  fbreaks,findices = freq_ticks(cohere.cort.aspect)

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
