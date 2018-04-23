function threshold_scales(cs,cutoff=2cycoct,buildup=1s)
    cs = cs[buildup .. last(times(cs))]
  threshold = mean(abs,cs)

  dims = (axisdim(cs,Axis{:scale}), axisdim(cs,Axis{:freq}))
  onestream = sum(max.(0,abs.(cs[:,0cycoct .. cutoff,:]) .-
                       threshold),dims)
  twostream = sum(max.(0,abs.(cs[:,cutoff + 1e-8cycoct .. last(scales(cs)),:]) .-
                       threshold),dims)

  onestreamax = AxisArray(squeeze(onestream,dims),Axis{:time}(times(cs)))
  twostreamax = AxisArray(squeeze(twostream,dims),Axis{:time}(times(cs)))

  onestreamax,twostreamax
end

function source_count_by_threshold(cs;cutoff=2cycoct,buildup=1s,window=1s,
                                   delta=0.25s)
  onestream,twostream = threshold_scales(cs,cutoff,buildup)
  onemean = map_window(mean,onestreamax,window,delta)
  twomean = map_window(mean,twostreamax,window,delta)

  AxisArray((twomean .> twomean) .+ 1,Axis{:time}(times(onemean)))
end
