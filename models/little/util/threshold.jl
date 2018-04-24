function source_count_by_threshold(cs;cutoff=2cycoct,buildup=1s,window=1s,
                                   delta=0.25s,return_intermediate=false)
  csa = AxisArray(abs.(cs),axes(cs)...)
  csa = csa[buildup .. last(times(cs))]
  threshold = mean(csa) + 2std(csa)

  dims = (axisdim(cs,Axis{:scale}), axisdim(cs,Axis{:freq}))
  onestream = sum(max.(0,csa[:,0cycoct .. cutoff,:] .- threshold),dims)
  twostream = sum(max.(0,csa[:,cutoff + 1e-8cycoct .. last(scales(cs)),:] .-
                       threshold),dims)

  onestream_ax = AxisArray(squeeze(onestream,dims),Axis{:time}(times(csa)))
  twostream_ax = AxisArray(squeeze(twostream,dims),Axis{:time}(times(csa)))

  onemean = map_window(mean,onestream_ax,window,delta)
  twomean = map_window(mean,twostream_ax,window,delta)

  onemean ./= maximum(onemean)
  twomean ./= maximum(twomean)

  result = AxisArray((onemean .> twomean) .+ 1,Axis{:time}(times(onemean)))

  if return_intermediate
    result, onemean, twomean
  else
    result
  end
end
