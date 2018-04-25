using AuditoryModel
using AxisArrays
const reasonable_response_maximum = 100

type ResponseOverflow <: Exception
  val::Float64
end

function bistable_scales(x,params)
  scales = cycoct.*2.0.^linspace(params[:scale_start],
                                 params[:scale_stop],
                                 params[:scale_N])

  noise_params = Dict(:τ_σ=>params[:τ_σ],:c_σ=>params[:c_σ])
  cortical_params = Dict(:scales=>scales)
  other_keys = [:τ_σ,:c_σ,:scale_start,:scale_stop,:scale_N,:W_m_σ]
  adapt_params = Dict(k => params[k] for k in setdiff(keys(params),other_keys))

  sp = audiospect(x,progressbar=false)
  cs = cortical(sp;cortical_params...,progressbar=false)

  scale_weights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                       axes(cs,Axis{:time}),
                            axes(cs,Axis{:scale}))

  weight_max = maximum(scale_weights)
  if weight_max > reasonable_response_maximum
    throw(ResponseOverflow(weight_max))
  end

  swn = drift(scale_weights;noise_params...,progressbar=false)
  swna,a,m = adaptmi(swn,W_m=scale_weighting(cs,params[:W_m_σ]),
                     shape_y = x -> max(0,x);
                     adapt_params...,progressbar=false)

  cs .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
end
