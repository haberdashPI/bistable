using VegaLite
using DataFrames

# export vraster

function vraster(z::AbstractMatrix;x=indices(z,1),y=indices(z,2),kwds...)
  Y = ones(x) .* y'
  X = x .* ones(y)'
  df = DataFrame(x = vec(X),y=vec(Y),z = vec(z))
  vraster(df;kwds...)
end

function vraster(df::DataFrame;value=:z,kwds...)
  if any(.!iszero.(imag.(df[value])))
    vraster_plot__complex(df,value=value;kwds...)
  else
    vraster_plot__real(df,value=value;kwds...)
  end
end

function vraster_plot__real(df::DataFrame;value=:z,x=:x,y=:y,maxbins=300,
                            limits=extrema(real.(df[value])),real_suffix=:_real)
  value_r = Symbol(string(value,real_suffix))
  df = copy(df)
  df[value_r] = real.(df[value])

  colorscheme = if any(df[value_r] .< 0)
    lim = maximum(abs.(limits))
    limits = (-lim,lim)
    "redblue"
  else
    "reds"
  end

  nx = length(unique(df[x]))
  ny = length(unique(df[y]))

  dataset(df) |>
    @vlplot(
      :rect,
      width=300, height=200,
      x={x, bin={maxbins=maxbins}},
      y={y, bin={maxbins=maxbins}},
      color="mean("*string(value)*")",
      config={range={ heatmap={ scheme=colorscheme } },
              view={ stroke="transparent" } }
     )
end

function vraster_plot__complex(df::DataFrame;value=:z,x=:x,y=:y,maxbins=300,
                              phase_suffix=:_phase,abs_suffix=:_abs)
  df = copy(df)
  z_phase = Symbol(string(value,phase_suffix))
  z_abs = Symbol(string(value,abs_suffix))
  df[z_phase] = angle.(df[value])
  df[z_abs] = abs.(df[value])

  nx = length(unique(df[x]))
  ny = length(unique(df[y]))

  dataset(df) |>
    @vlplot(
      :rect,
      width=300, height=200,
      x={x, bin={maxbins=maxbins}},
      y={y, bin={maxbins=maxbins}},
      color="mean("*string(z_phase)*")",
      opacity="mean("*string(z_abs)*")",
      config={range={ heatmap={ scheme="sinebow" } },
              view={ stroke="transparent" } }
     )
end
