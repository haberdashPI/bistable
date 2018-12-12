using Gadfly
using CSV
using DependentBootstrap
using DataFramesMeta
using ShiftedArrays
using PlotAxes
PlotAxes.set_backend!(:gadfly)

function packing(x;maxpad=true)
    vals = sort!(unique(x))
    pos = [0; cumsum([1.5; fill(1,length(vals)-3); maxpad ? 1.5 : 1])]
    vals, pos
end

function packaxes(x;kwds...)
    p = Dict(zip(packing(x;kwds...)...))
    [p[xi] for xi in x]
end

function packaxes_invfn(x;kwds...)
   p,v = packing(x;kwds...)
   p = Dict(zip(v,p))
   xi -> p[xi]
end

function rename_levels_for(df,vals;suffixes=[:c_a,:c_m])
  for suf in suffixes
    df[suf] = zero(df[Symbol("f_"*string(suf))])
  end
  df[:level] = "unknown"
  for (prefix,level) in [("f_","Peripheral"),("s_","Cortical"),("t_","Object")]
    rows = df[Symbol(prefix*"c_σ")] .> 0
    df[rows,:level] = level
    for suffix in suffixes
      df[rows,suffix] = df[rows,Symbol(prefix*string(suffix))]
    end
  end
  df[[suffixes;:level;vals]]
end

DataFramesMeta.linq(::DataFramesMeta.SymbolParameter{:rename_levels}, df, vals) = :(rename_levels($df,$vals))

function colorscale(smap;reverse=false,minvalue,maxvalue,colorstart=minvalue,
                    colorstop=maxvalue,
                    midvalue=(colorstop-colorstart)/2+colorstart)
    scale(x,min,max) = (x - min)/(max - min)
    scalei(x,min,max) = x*(max - min) + min
    fi(x) = clamp(scalei(x,minvalue,maxvalue),minvalue,maxvalue)
    c(x) = scale(x,colorstart,colorstop)
    colors = !reverse ? colormap(smap,mid=c(midvalue)) :
      reverse!(colormap(smap,mid=1-c(midvalue)))

    function(x)
        colors[1+clamp(floor(Int,99*c(fi(x))),0,99)]
    end
end

function select_mask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                     stop_time=20s,kwds...)
  selections = select_params(params;Δf=Δf,kwds...)
  if length(selections) > 1
    error("Ambiguous parameter selection, mutliple matches.")
  end
  selection = selections[1]
  masks = []
  for_results_in(joinpath(datadir,"data")) do entry
    if entry["pindex"] == selection[1]
      push!(masks,entry["mask"])
    end
  end
  audiospect(masks[simulation],settings)[start_time .. stop_time]
end

function unzip(xs)
  fns = map(enumerate(first(xs))) do (i,el)
    ys -> map(y -> y[i],ys)
  end |> Tuple
  map(fns) do fn
    fn(xs)
  end
end

# okay, the below is *clearly* wrong somehow, and I don't
# know how to fix it right now; abandon for the moment

# """
#     qqci(x,y,indices)

# Compute a confidence interval around each quantile of x
# over the samples of y selected by indices. Indicies
# is an iterable where each element is an array of indices
# to select: e.g. the output of dbootinds.
# """
# function qqci(x,y,indices;alpha=0.05)
#   x, yvals = qqdistboot(x,y,indices)
#   ylower, yupper = mapslices(yvals,dims=2) do samples
#     quantile(samples,(alpha/2,1-alpha/2))
#   end |> unzip

#   x,vec(ylower),vec(yupper)
# end

# function qqdistboot(x,y,indices)
#   nx = length(x)
#   yvals = Array{float(eltype(x))}(undef,nx,length(indices))

#   for (i,ixs) in enumerate(indices)
#     qx, qy, order = qqdata(x,y[ixs])
#     lastxi = 1
#     for (j,xi) in enumerate(order)
#       if xi <= nx
#         yvals[xi,i] = qy[j]
#         lastxi = xi
#       else
#         yvals[lastxi,i] = qy[j]
#       end
#     end
#   end

#   cumsum(sort(x)) ./ sum(x), yvals
# end

function qqdata(x,y)
  nx = length(x)
  ny = length(y)
  sumx = sum(x)
  sumy = sum(y)
  order = sortperm([x; y])
  qx,qy = mapreduce((x,a) -> push!(x,x[end].+a),order,init=[(0.0,0.0)]) do i
    i <= nx ? (x[i]/sumx,0.0) : (0.0,y[i-nx]/sumy)
  end |> unzip

  qx, qy, order
end

function plot_fit(df,params;separate_plots=false,showci=true,kwds...)
  df,params = select_data(df,params;kwds...)
  hdata = human_data()

  sim = by(stream_summary(df,params),:st) do g
    DataFrame([findci(g.streaming)])
  end
  sim[:experiment] = "simulation"

  human = by(hdata.stream,:st) do g
    DataFrame([findci(g.streaming)])
  end
  human[:experiment] = "human"

  stream = vcat(sim,human)
  stream = @transform(stream,
                      pos = @.(log2(:st) -
                               ifelse(:experiment == "human",0.05,-0.05)));

  splot = plot(stream,x=:pos,y=:mean,ymin=:lowerc,ymax=:upperc,
     color=:experiment,shape=:experiment,
     Guide.shapekey(pos=[log2(3),1.4]),
     Guide.xticks(ticks=log2.([3,6,12])),
     Scale.color_discrete_manual("lightgray","darkgray"),
     Scale.x_continuous(labels=x -> string(floor(Int,2.0^x))),
     Scale.y_continuous(labels=x -> string(floor(Int,100*x))),
     Guide.xlabel("Δf (semitones)",orientation=:horizontal),
     Guide.ylabel("% streaming",orientation=:vertical),
     Coord.cartesian(xmin=log2(3)-0.25,xmax=log2(12)+0.25),
     Geom.point,Geom.line,Geom.errorbar,
     Theme(discrete_highlight_color=x->"black"))

  hlens = normlength(hdata.lengths)
  slens = normlength(length_summary(df,params))

  qh,qs = qqdata(hlens,slens)
  nsamples = div(length(hlens),N_for_pressnitzer_hupe_2006)
  mean = DataFrame(qhuman = qh, qsim = qs)

  # cisim,cihuman_lower,cihuman_upper =
  #   qqci(slens,hlens, dbootinds(hlens,numobsperresample=nsamples,
  #                               numresample=1000))
  # bounds = DataFrame(human_lower = cihuman_lower,human_upper = cihuman_upper,
  #                    sim = cisim)

  hplot = plot(mean,x=:qhuman,y=:qsim,xintercept=[0],slope=[1],
               Geom.line,Geom.abline,
               Theme(default_color="grey"),Coord.cartesian(xmax=1,ymax=1))

  # hplot = plot(layer(mean,x=:qhuman,y=:qsim,xintercept=[0],slope=[1],
  #                    Geom.line,Geom.abline,
  #                    Theme(default_color="grey")),
  #              layer(bounds,x=:human_lower,y=:sim,Geom.line,
  #                    Theme(default_color="lightgrey")),
  #              layer(bounds,x=:human_upper,y=:sim,Geom.line,
  #                    Theme(default_color="lightgrey")))

#   lens = [DataFrame(nlength=hlens,experiment="human");
#           DataFrame(nlength=slens,experiment="simulation")]
#   hplot = plot(lens,x=:nlength,color=:experiment,
#                Guide.colorkey(pos=[0.65*Gadfly.w,-0.3*Gadfly.h]),
#                Scale.color_discrete_manual("darkgray","lightgray"),
#                Coord.cartesian(xmin=0,xmax=20),
#                Guide.xlabel("log-z-scored length",
#                             orientation=:horizontal),
#                Guide.ylabel("density",orientation=:vertical),
#                Geom.histogram(position=:dodge,bincount=60,density=true),
#                Theme(discrete_highlight_color=x->"black",
#                      bar_highlight=x->"black"))

  if separate_plots
    splot,hplot
  else
    hstack(splot,hplot)
  end
end

function plot_lengths((len,value))
  dfl = DataFrame(value=value,time=cumsum(len));
  dfl = @transform(dfl,lagtime = lag(:time,default=0.0),ymin=-1,ymax=1);

  plot(dfl,xmax=:time,ymin=:ymin,xmin=:lagtime,ymax=:ymax,color=:value,
       Scale.color_discrete_manual("black","lightgray"),
       Geom.rect)
end

function plot_fitmask(data,params,settings;Δf=6,simulation=1,start_time=0s,
                      stop_time=20s,kwds...)
  fit = plot_fit(data,params;kwds...)
  mask = plot_mask(data,params,settings;Δf=Δf,simulation=simulation,
                   start_time=start_time,stop_time=stop_time,kwds...)
  vstack(fit,mask)
end

function plot_mask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                   stop_time=20s,kwds...)
  selection = select_params(params;Δf=Δf,kwds...)
  mask = select_mask(df,params,settings;Δf=Δf,simulation=simulation,
                     start_time=start_time,stop_time=stop_time,kwds...)
  input = audiospect_stimulus(params[selection,:],settings)
  input = input[start_time .. stop_time]
  band = plot_lengths(percept_lengths(mask,input,settings))

  spect=plot(asplotable(mask,quantize=(200,128))[1],
             x=:time,y=:logfreq,color=:value,Geom.rectbin,
             Coord.cartesian(xmin=ustrip(start_time),xmax=ustrip(stop_time)),
             Scale.color_continuous(colormap=Scale.lab_gradient("white","black")))

  vstack(band,spect)
end
