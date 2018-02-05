using RCall
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

  sindices = 1:2:length(scales(cohere.cort))
  sbreaks = round.(scales(cohere.cort)[sindices],2)
  fbreaks,findices = freq_ticks(cohere.cort.aspect,u[:,:,1])

  p = raster_plot(df,value=:response,x=:scale_index,y=:freq_bin)

R"""
    scale_color_brewer(palette='Set1',name=$(string(label[1]))) +
    coord_cartesian(ylim=c(1.1,0)) + ylab('lambda / var(x)') +
    xlab('time (s)')
"""
end

  
