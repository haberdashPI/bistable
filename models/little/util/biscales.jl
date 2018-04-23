using AuditoryModel
using AxisArrays
const reasonable_response_maximum = 100

type ResponseOverflow <: Exception
  val::Float64
end

function bistable_scales(x,params)
  noise_params = Dict(:τ_σ=>params[:τ_σ],:c_σ=>params[:c_σ])

  adapt_params = Dict(:c_m=>params[:c_m],
                      :τ_m=>params[:τ_m],
                      :c_e=>params[:c_e],
                      :τ_e=>params[:τ_e],
                      :c_a=>params[:c_a],
                      :τ_a=>params[:τ_a],

                      :c_e=>params[:c_e],
                      :τ_e=>params[:τ_e],
                      :c_d=>params[:c_d],
                      :τ_d=>params[:τ_d],

                      :α=>params[:α])

  cortical_params = Dict(:scales=>params[:scales])

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
  swna,a,m = adaptmi(swn,W_m=scale_weighting2(cs,params[:W_m_σ]),
                     shape_y = x -> max(0,x);
                     adapt_params...,progressbar=false)

  cs .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
end
