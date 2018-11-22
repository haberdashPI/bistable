using Gadfly
using CSV
using DependentBootstrap
using DataFramesMeta
using ShiftedArrays
using PlotAxes

function select_params(params;kwds...)
  condition = trues(size(params,1))
  for (var,val) in pairs(kwds)
    condition .&= abs.(params[var] .- val) .<= 0.1
  end
  params[condition,:pindex]
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


function stream_dfs(df,params;keep_simulations=false,findci=true,
                    bound=true,bound_threshold=0.8,kwds...)
  selection = select_params(params;kwds...)
  if length(selection) != 3
    error("Expected three parameter entires, one for each Δf.",
          "\nInstead found the entires: ",string(selection),
          "\nKeyword selection: ",string(kwds))
  end
  dfstream_ind = @linq df |>
    where(in.(:pindex,Ref(selection))) |>
    by([:pindex,:created],
       st = params[first(:pindex),:].Δf,
       streaming = streamprop(:percepts,:length,bound=bound,
                              threshold=bound_threshold),
       streaming_unbound = streamprop(:percepts,:length,bound=false))

  dfstream = by(dfstream_ind,:st) do dfind
    if any(!ismissing,dfind[:streaming])
      @with dfind begin
        streaming = collect(skipmissing(:streaming))
        if findci
          if length(streaming) > 3 &&
              any(x != streaming[1] for x in streaming)
            cint = dbootconf(streaming)
            DataFrame(mean = mean(streaming),
                      lowerc = cint[1],
                      upperc = cint[2])
          else
            DataFrame(mean = mean(streaming),
                      lowerc = minimum(streaming),
                      upperc = maximum(streaming))
          end
        else
          DataFrame(mean = mean(streaming))
        end
      end
    else
      # if we don't have any (or very little) bistability at all
      # the analysis above will fail: we still want to report the
      # dominant percept (but no estimate of bounds)
      val = @with(dfind,mean(:streaming_unbound))
      if findci
        DataFrame(mean = val, lowerc=val, upperc=val)
      else
        DataFrame(mean = val)
      end
    end
  end

  dfstream.experiment = "simulation"

  if keep_simulations
    dfstream, dfstream_ind
  else
    dfstream
  end
end

function stream_df(df,params;findci=true,bound=true,bound_threshold=0.8,kwds...)
  dfs = stream_dfs(df,params;findci=findci,bound=bound,
                   bound_threshold=bound_threshold,kwds...)
  dfh = stream_dfh(findci=findci)
  delete!(dfh,:rms)
  vcat(dfs,dfh)
end

function length_dfh()
  ph = CSV.read(joinpath("..","data","pressnitzer_hupe",
                         "pressnitzer_hupe_inferred.csv"))
  DataFrame(length = ph[:length],experiment="human",
            nlength = exp.(zscore(log.(ph[:length]))));
end

function length_dfs(df,params;kwds...)
  selection = select_params(params;Δf=6,kwds...)
  dfs = @where(df,(:pindex .== selection))
  DataFrame(length = dfs[:length],experiment="simulation",
            nlength = exp.(zscore(log.(dfs[:length]))));
end

length_df(df,params;kwds...) = vcat(length_dfh(), length_dfs(df,params;kwds...))

function plot_fitmask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                      stop_time=20s,kwds...)
  fit = plot_fit(df,params;kwds...)
  mask = plot_mask(df,params,settings;Δf=Δf,simulation=simulation,
                   start_time=start_time,stop_time=stop_time,kwds...)
  vstack(fit,mask)
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

function handlebound(fn,seconds;bound=true,threshold=0.8)
    if bound && length(seconds) < 3
        return missing
    end

    if !bound || (sum(seconds[2:end-1]) > threshold*sum(seconds))
        fn(1:length(seconds))
    else
        fn(2:length(seconds)-1)
    end
end

function streamprop(percepts,seconds;kwds...)
    handlebound(seconds;kwds...) do range
        sum(seconds[range][percepts[range] .== 2]) / sum(seconds[range])
    end
end

function stim_per_second(seconds;kwds...)
    handlebound(seconds;kwds...) do range
        length(range) / sum(seconds[range])
    end
end
