using Gadfly
using CSV
using DependentBootstrap
using DataFramesMeta
using ShiftedArrays
using PlotAxes

function packing(x)
    vals = sort!(unique(x))
    pos = [0; cumsum([1.5; fill(1,length(vals)-3); 1.5])]
    vals, pos
end

function packaxes(x)
    p = Dict(zip(packing(x)...))
    [p[xi] for xi in x]
end

function packaxes_invfn(x)
   p,v = packing(x)
   p = Dict(zip(v,p))
   xi -> p[xi]
end

function rename_levels_for(df,vals)
    df[:c_m] = NaN
    df[:c_a] = NaN
    df[:level] = "unknown"
    @byrow! df begin
        if :f_c_σ > 0
           :level = "Peripheral"
            :c_m = :f_c_m
            :c_a = :f_c_a
        elseif :s_c_σ > 0
            :level = "Cortical"
            :c_m = :s_c_m
            :c_a = :s_c_a
        elseif :t_c_σ > 0
            :level = "Object"
            :c_m = :t_c_m
            :c_a = :t_c_a
        end
    end
    df[[:c_m;:c_a;:level;vals]]
end

DataFramesMeta.linq(::DataFramesMeta.SymbolParameter{:rename_levels}, df, vals) = :(rename_levels($df,$vals))

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

function plot_lengths((len,value))
  dfl = DataFrame(value=value,time=cumsum(len));
  dfl = @transform(dfl,lagtime = lag(:time,default=0.0),ymin=-1,ymax=1);

  plot(dfl,xmax=:time,ymin=:ymin,xmin=:lagtime,ymax=:ymax,color=:value,
       Scale.color_discrete_manual("black","lightgray"),
       Geom.rect)
end

function plot_mask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                   stop_time=20s,kwds...)
  selection = select_params(params;Δf=Δf,kwds...)
  mask = select_mask(df,params,settings;Δf=Δf,simulation=simulation,
                     start_time=start_time,stop_time=stop_time,kwds...)
  input = audiospect_stimulus(params[selection,:],settings)
  input = input[start_time .. stop_time]
  plot_lengths(percept_lengths(mask,input,settings))

  spect=plot(asplotable(mask,quantize_size=(200,128))[1],
             x=:time,y=:logfreq,color=:value,Geom.rectbin,
             Coord.cartesian(xmin=ustrip(start_time),xmax=ustrip(stop_time)),
             Scale.color_continuous(colormap=Scale.lab_gradient("white","black")))

  vstack(band,spect)
end
