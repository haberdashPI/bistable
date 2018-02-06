using RCall
export rplot, collapsed_scale_plot

rplot(cort::CorticalModel,y::AbstractVector) = rplot(cort,cort(y))
rplot(cort::CorticalModel,y::AbstractMatrix) = rplot(cort,cort(y))

function nearin(xin,xs)
  _,inds = findmin(abs.(vec(xin) .- vec(xs)'),2)
  cols = map(ii -> ii[2],CartesianRange((length(xin),length(xs))))
  if length(inds) > 0
    cols[inds[:,1]]
  else
    Array{Int}(0)
  end
end

function findnear(x,nearby)
  indices = sort(unique(nearin(filter(!isnan,nearby),
                                filter(!isnan,x))))
  if any(isnan(nearby))
    [NaN; x[indices]], [find(isnan,x); indices]
  else
    x[indices], indices
  end
end

function rplot(cort::CorticalModel,y;rates=cort.rates,scales=cort.scales)
  rates, rindices = findnear(cort.rates,rates)
  scales, sindices = findnear(cort.scales,scales)

  if rates != cort.rates || scales != cort.scales
    @show rindices
    @show sindices
  end

  y = y[:,rindices,sindices,:]

  ixs = CartesianRange(size(y))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(y),
                 time = ustrip(vec(times(cort,y)[at(ixs,1)])),
                 rate = vec(rates[at(ixs,2)]),
                 scale = vec(scales[at(ixs,3)]),
                 freq_bin = vec(at(ixs,4)))

  fbreaks,findices = freq_ticks(cort.aspect)
  p = raster_plot(df,value=:response,x=:time,y=:freq_bin)

R"""

  library(ggplot2)

  scalestr = function(x){
    ifelse(!is.nan(x),sprintf("Scale: %3.2f cyc/oct",x),"All Scales")
  }
  ratestr = function(x){
    ifelse(!is.nan(x),sprintf("Rate: %5.2f Hz",x),"All Rates")
  }

  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr($(sort(scales))))
  }
  ordered_rates = function(x){
    factor(ratestr(x),levels=ratestr($(sort(rates))))
  }

  $p +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(ordered_scales(scale) ~ ordered_rates(rate))

"""
end

# function collapsed_scale_plot(cort,m,range=nothing)
#   sm = m[:,1,:,1]
#   iis = collect(CartesianRange(size(sm)))
#   df = DataFrame(level = vec(sm),
#                  time = vec(map(x -> ustrip(times(cort,m)[x[1]]),iis)),
#                  scale = vec(map(x -> ustrip(scales(cort)[x[2]]),iis)))

#   @show unique(df[:scale])
# R"""
#   library(RColorBrewer)
#   pal = brewer.pal(5,'Reds')[2:5]
#   p = ggplot($df,aes(x=time,y=level,group=scale,
#                  color=factor(round(scale,1)),linetype=factor(round(scale,1)))) +
#     geom_line(color='black',linetype='solid',size=1.2) + geom_line() +
#     scale_color_manual(values = rep(pal,each=3)[1:11],name="Scale") +
#     scale_linetype_manual(values =
#       rep(c("dotdash","longdash","solid"),4)[1:11],name="Scale")
# """
#   if range != nothing
# R"""
#       p = p + coord_cartesian(ylim=c($(first(range)),$(last(range))))
# """
#   end

#   R"p"
# end

function collapsed_scale_plot(cort,data;name="response",range=nothing)
  data = data[:,1,:,1]
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(cort,data)[at(ixs,1)])),
                 scale_bin = vec(at(ixs,2)))

  sbreaks = 1:2:length(scales(cort))
  slabs = string.(round.(scales(cort)[sbreaks],2))

  p = raster_plot(df,value=:response,x=:time,y=:scale_bin,name=name)

R"""

  library(ggplot2)

  $p +
    scale_y_continuous(labels=$slabs,breaks=$sbreaks) +
    ylab('Scale (cyc/oct)') + xlab('Time (s)')

"""
end


function plot_scales2(cort,data::Array{<:Complex};name="response",range=nothing)
  data = data[:,1,:,1]
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(cort,data)[at(ixs,1)])),
                 scale_bin = vec(at(ixs,2)))

  sbreaks = 1:2:length(scales(cort))
  slabs = string.(round(scales(cort)[sbreaks],2))

  p = raster_plot(df,value=:response,x=:time,y=:scale_bin)

R"""

  library(ggplot2)

  $p +
    scale_y_continuous(labels=$slabs,breaks=$sbreaks) +
    ylab('Scale (cyc/oct)') + xlab('Time (s)')

"""
end
