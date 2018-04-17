using AuditoryModel

function bistable_scales(x,noise_params,adapt_params,cortical_params)
  sp = audiospect(x,progressbar=false)
  cs = cortical(sp;cortical_params...,progressbar=false)

  scale_weights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                       axes(cs,Axis{:time}),
                       axes(cs,Axis{:scale}))

  swn = drift(scale_weights;progressbar=false,noise_params...)
  swna,a,m = adaptmi(swn,W_m=scale_weighting2(cs,1.0),shape_y = x -> max(0,x);
                     progressbar=false,adapt_params...)

  cs .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
end
