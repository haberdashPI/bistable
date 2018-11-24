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

function plot_fit(df,params;kwds...)
  lens = length_df(df,params;kwds...)

  stream = @transform(stream_df(df,params;kwds...),
                      pos = @.(log(:st)/log(2) -
                               ifelse(:experiment == "human",0.05,-0.05)));
  splot = plot(stream,x=:pos,y=:mean,ymin=:lowerc,ymax=:upperc,
     color=:experiment,shape=:experiment,
     Geom.point,Geom.line,Geom.errorbar,Coord.cartesian(xmin=1,xmax=4))

  hplot = plot(lens,x=:nlength,color=:experiment,
               Coord.cartesian(xmin=0,xmax=10),
               Geom.histogram(position=:dodge,bincount=60,density=true))

  hstack(splot,hplot)
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
       Geom.rect)

  spect=plot(asplotable(mask,quantize_size=(200,128))[1],
       x=:time,y=:logfreq,color=:value,Geom.rectbin,
       Coord.cartesian(xmin=ustrip(start_time),xmax=ustrip(stop_time)),
       Scale.color_continuous(colormap=Scale.lab_gradient("white","red")))

  vstack(band,spect)
end
