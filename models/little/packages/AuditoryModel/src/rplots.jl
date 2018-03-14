using RCall
using DataFrames
using JLD

export rplot, collapsed_scale_plot
import Colors: RGB

const cmap = Dict{Symbol,Vector{RGB}}()
jldopen(joinpath(@__DIR__,"..","data","colormaps.jld")) do file
  cmap[:D1] = read(file,"D1")
  cmap[:reds] = read(file,"reds")
  cmap[:C6] = read(file,"C6")
end

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
    "#".*hex.(cmap[:D1])
  else
    "#".*hex.(cmap[:reds])
  end

R"""

  ggplot($df,aes_string(x=$(string(x)),y=$(string(y)),fill=$(string(value_r)))) +
    geom_raster() +
    scale_fill_gradientn(colors=$colormap,limits=$(collect(limits)),name="x")

"""
end
function raster_plot__complex(df::DataFrame;value=:z,x=:x,y=:y,
                              phase_suffix=:_phase,abs_suffix=:_abs)
  colormap = "#".*hex.(cmap[:C6])
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

function rplot(as::AuditorySpectrogram)
  ixs = CartesianRange(size(as))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(as),
                 time = vec(ustrip(times(as)[at(ixs,1)])),
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

function rplot(cort::CorticalRates;rates=AuditoryModel.rates(cort))
  cort = cort[:,atvalue.(rates),:]
  @show size(cort)
  @show typeof(cort)
  ixs = CartesianRange(size(cort))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(cort),
                 time = ustrip.(vec(times(cort)[at(ixs,1)])),
                 rate = ustrip.(vec(rates[at(ixs,2)])),
                 freq_bin = vec(at(ixs,3)))

  fbreaks,findices = freq_ticks(cort)
  p = raster_plot(df,value=:response,x=:time,y=:freq_bin)

R"""

  library(ggplot2)

  ratevals = $(ustrip.(rates))
  ratestr = function(x){
    ifelse(!is.nan(x),sprintf("Rate: %5.2f Hz",x),"All Rates")
  }

  ordered_rates = function(x){
    factor(ratestr(x),levels=ratestr(ratevals))
  }

  $p +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(~ordered_rates(rate))

"""
end

function rplot(cort::CorticalScales;scales=AuditoryModel.scales(cort))
  cort = cort[:,atvalue.(scales),:]
  ixs = CartesianRange(size(cort))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(cort),
                 time = ustrip.(vec(times(cort)[at(ixs,1)])),
                 scale = ustrip.(vec(scales[at(ixs,2)])),
                 freq_bin = vec(at(ixs,3)))

  fbreaks,findices = freq_ticks(cort)
  p = raster_plot(df,value=:response,x=:time,y=:freq_bin)

R"""

  library(ggplot2)
  scalevals = $(ustrip.(scales))
  scalestr = function(x){
    ifelse(!is.nan(x),sprintf("Scale: %3.2f cyc/oct",x),"All Scales")
  }

  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr(scalevals))
  }

  $p +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(ordered_scales(scale) ~ .)

"""
end


function rplot(cort::Cortical;rates=AuditoryModel.rates(cort),
               scales=AuditoryModel.scales(cort))
  cort = cort[:,atvalue.(rates),atvalue.(scales),:]
  ixs = CartesianRange(size(cort))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(cort),
                 time = ustrip.(vec(times(cort)[at(ixs,1)])),
                 rate = ustrip.(vec(rates[at(ixs,2)])),
                 scale = ustrip.(vec(scales[at(ixs,3)])),
                 freq_bin = vec(at(ixs,4)))

  fbreaks,findices = freq_ticks(cort)
  p = raster_plot(df,value=:response,x=:time,y=:freq_bin)

R"""

  library(ggplot2)

  scalevals = $(ustrip.(scales))
  ratevals = $(ustrip.(rates))
  scalestr = function(x){
    ifelse(!is.nan(x),sprintf("Scale: %3.2f cyc/oct",x),"All Scales")
  }
  ratestr = function(x){
    ifelse(!is.nan(x),sprintf("Rate: %5.2f Hz",x),"All Rates")
  }

  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr(scalevals))
  }
  ordered_rates = function(x){
    factor(ratestr(x),levels=ratestr(ratevals))
  }

  $p +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (kHz)') + xlab('Time (s)') +
    facet_grid(ordered_scales(scale) ~ ordered_rates(rate))

"""
end

function collapsed_scale_plot(cort;range=nothing)
  scaledim = axisdim(AxisArray(cort),Axis{:scale})
  timedim = axisdim(AxisArray(cort),Axis{:time})

  dims = Tuple(find(indexin(1:ndims(cort),[scaledim,timedim]) .== 0))
  y = squeeze(mean(abs.(AxisArray(cort).data),dims),dims)
  ixs = CartesianRange(size(y))
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(y),
                 time = vec(ustrip(times(cort)[at(ixs,1)])),
                 scale_bin = vec(at(ixs,2)))

  sbreaks = 1:2:length(scales(cort))
  slabs = string.(round.(ustrip.(uconvert.(cycoct,scales(cort)))[sbreaks],2))

  p = raster_plot(df,value=:response,x=:time,y=:scale_bin)

R"""

  library(ggplot2)

  $p +
    scale_y_continuous(labels=$slabs,breaks=$sbreaks) +
    ylab('Scale (cyc/oct)') + xlab('Time (s)')

"""
end
