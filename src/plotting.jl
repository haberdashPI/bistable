using Gadfly
using CSV
using DependentBootstrap
using DataFramesMeta
using ShiftedArrays
using PlotAxes

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

function plot_fit(df,params;separate_plots=false,kwds...)
  lens = @where(length_df(df,params;kwds...),
                :nlength .<= maximum(normalized_hist_range))

  stream = @transform(stream_df(df,params;kwds...),
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

  hplot = plot(lens,x=:nlength,color=:experiment,
               Guide.colorkey(pos=[0.65*Gadfly.w,-0.3*Gadfly.h]),
               Scale.color_discrete_manual("darkgray","lightgray"),
               Coord.cartesian(xmin=0,xmax=20),
               Guide.xlabel("log-z-scored length",
                            orientation=:horizontal),
               Guide.ylabel("density",orientation=:vertical),
               Geom.histogram(position=:dodge,bincount=30,density=true),
               Theme(discrete_highlight_color=x->"black",
                     bar_highlight=x->"black"))

  if separate_plots
    splot,hplot
  else
    hstack(splot,hplot)
  end
end

function plot_mask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                   stop_time=20s,kwds...)
  selection = select_params(params;Δf=Δf,kwds...)
  mask = select_mask(df,params,settings;Δf=Δf,simulation=simulation,
                     start_time=start_time,stop_time=stop_time,kwds...)
  input = audiospect_stimulus(params[selection,:],settings)
  input = input[start_time .. stop_time]
  l,v = percept_lengths(mask,input,settings)
  dfl = DataFrame(value=v,time=ustrip(start_time).+cumsum(l));
  dfl = @transform(dfl,lagtime = lag(:time,default=float(ustrip(start_time))),
                   ymin=-1,ymax=1);

  band=plot(dfl,xmax=:time,ymin=:ymin,xmin=:lagtime,ymax=:ymax,color=:value,
               Scale.color_discrete_manual("black","lightgray"),
               Geom.rect)

  spect=plot(asplotable(mask,quantize_size=(200,128))[1],
       x=:time,y=:logfreq,color=:value,Geom.rectbin,
       Coord.cartesian(xmin=ustrip(start_time),xmax=ustrip(stop_time)),
       Scale.color_continuous(colormap=Scale.lab_gradient("white","black")))

  vstack(band,spect)
end
