using AuditoryModel
using AxisArrays
const reasonable_response_maximum = 100

type ResponseOverflow <: Exception
  val::Float64
end

function meanabs(A,n)
  result = reducedim((x,y) -> x+abs(y),A,n,zero(real(eltype(A))))
  result ./= size(A,n)
  squeeze(result,n)
end

function bistable_scales(cs,params)
  noise_params = Dict(:τ_σ=>params[:τ_σ],:c_σ=>params[:c_σ])
  other_keys = [:τ_σ,:c_σ,:W_m_σ]
  adapt_params = Dict(k => params[k] for k in setdiff(keys(params),other_keys))

  scale_weights = AxisArray(meanabs(cs,axisdim(cs,Axis{:freq})),
                            axes(cs,Axis{:time}), axes(cs,Axis{:scale}))

  swn = drift(scale_weights;noise_params...,progressbar=false)
  swna,a,m = adaptmi(swn;W_m=scale_weighting(cs,params[:W_m_σ]),
                     shape_y = x -> max(0,x),progressbar=false,
                     adapt_params...)

  weight_max = maximum(swna)
  if weight_max > reasonable_response_maximum
    throw(ResponseOverflow(weight_max))
  end

  cs .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs).*im)
end
