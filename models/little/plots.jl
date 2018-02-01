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

function raster_plot__real(df::DataFrame;value=:z,x=:x,y=:y,name=string(value),
                           limits=extrema(df[value]))

  colormap = if any(df[value] .< 0)
    lim = maximum(abs.(limits))
    limits = (-lim,lim)
    "#".*hex.(RGB.(cmap("D1")))
  else
    "#".*hex.(RGB.(Colors.colormap("Reds")))
  end

R"""

  ggplot($df,aes_string(x=$(string(x)),y=$(string(y)),fill=$(string(value)))) +
    geom_raster() +
    scale_fill_gradientn(colors=$colormap,limits=$(collect(limits)),name=$name)

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
