using RCall
export rplot, collapsed_scale_plot
import PerceptualColourMaps: cmap
import Colors: RGB

R"library(ggplot2)"

function raster_plot(z::Matrix;x=indices(z,1),y=indices(z,2))
  Y = ones(x) .* y'
  X = x .* ones(y)'
  df = DataFrame(x = vec(X),y=vec(Y),z = vec(z))
  raster_plot(df)
end

function raster_plot(df::DataFrame;value=:z,kwds...)
  if any(.!iszero.(imag.(df[value])))
    raster_plot__complex(df,value=value;kwds...)
  else
    raster_plot__real(df,value=value;kwds...)
  end
end

function raster_plot__real(df::DataFrame;value=:z,x=:x,y=:y,
                           limits=extrema(real.(df[value])),real_suffix=:_real)
  value_r = Symbol(string(value,real_suffix))
  df = copy(df)
  df[value_r] = real.(df[value])

  colormap = if any(df[value_r] .< 0)
    lim = maximum(abs.(limits))
    limits = (-lim,lim)
    "#".*hex.(RGB.(cmap("D1")))
  else
    "#".*hex.(RGB.(Colors.colormap("Reds")))
  end

R"""

  ggplot($df,aes_string(x=$(string(x)),y=$(string(y)),fill=$(string(value_r)))) +
    geom_raster() +
    scale_fill_gradientn(colors=$colormap,limits=$(collect(limits)),name="x")

"""
end
function raster_plot__complex(df::DataFrame;value=:z,x=:x,y=:y,
                              phase_suffix=:_phase,abs_suffix=:_abs)
  colormap = "#".*hex.(RGB.(cmap("C6")))
  df = copy(df)
  z_phase = Symbol(string(value,phase_suffix))
  z_abs = Symbol(string(value,abs_suffix))
  df[z_phase] = angle.(df[value])
  df[z_abs] = abs.(df[value])

R"""

  ggplot($df,aes_string(x=$(string(x)),y=$(string(y)),
                        fill=$(string(z_phase)),alpha=$(string(z_abs)))) +
    geom_raster() +
    scale_fill_gradientn(colors=$colormap,limits=c(-pi-0.01,pi+0.01),
                         breaks=c(-pi,0,pi),
                         labels=c(expression(-pi),expression(0),
                                  expression(+pi)),
                         name = expression(Arg(x)))+
    scale_alpha_continuous(range=c(0,1),name="|x|")

"""
end

function aplot(data) =
  @assert haskey(data,"kind") "Missing 'kind' metadata"
  aplot(data,data["kind"])
end

aplot(data,as::AuditorySpectrogram) = rplot(as,audiospect(data,as))
function aplot(data::ImageMetaAxis,as::AuditorySpectrogram)
  @assert haskey(data,"kind") "Missing 'kind' metadata"
  @assert as == data["kind"] "Metadata must match spectrogram parameters"

  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(as,data)[at(ixs,1)])),
                 freq_bin = vec(at(ixs,2)))
  fbreaks,findices = freq_ticks(as)
  p = raster_plot(df,value=:response,x=:time,y=:freq_bin)
R"""

  library(ggplot2)

  $p +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (Hz)') + xlab('Time (s)')

"""
end

rplot(cort::CorticalModel,y::AbstractVector;kwds...) = rplot(cort,cort(y);kwds...)
rplot(cort::CorticalModel,y::AbstractMatrix;kwds...) = rplot(cort,cort(y);kwds...)

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
  if any(isnan.(nearby))
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

function collapsed_scale_plot(cort,data;range=nothing)
  data = data[:,1,:,1]
  ixs = CartesianRange(size(data))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(data),
                 time = vec(ustrip(times(cort,data)[at(ixs,1)])),
                 scale_bin = vec(at(ixs,2)))

  sbreaks = 1:2:length(scales(cort))
  slabs = string.(round.(scales(cort)[sbreaks],2))

  p = raster_plot(df,value=:response,x=:time,y=:scale_bin)

R"""

  library(ggplot2)

  $p +
    scale_y_continuous(labels=$slabs,breaks=$sbreaks) +
    ylab('Scale (cyc/oct)') + xlab('Time (s)')

"""
end


function plot_scales2(cort,data::AbstractArray{<:Complex};name="response",
                      range=nothing)
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
