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

function bistable_scales(cs,params;intermediate_results=false,
                        progressbar=false)
  if params[:condition] != :scales
    return cs
  end

  noise_params = Dict(:τ_σ=>params[:τ_σ],:c_σ=>params[:c_σ])
  adapt_params = Dict(
    :c_m => params[:c_m], :τ_m => params[:τ_m],
    :c_a => params[:c_a], :τ_a => params[:τ_a],
    :τ_x => params[:τ_x], :c_x => params[:c_x],
    :c_n => params[:c_n], :τ_n => params[:τ_n]
  )

  scale_weights = AxisArray(meanabs(cs,axisdim(cs,Axis{:freq})),
                            axes(cs,Axis{:time}), axes(cs,Axis{:scale}))

  swn = drift(scale_weights;noise_params...,progressbar=progressbar)
  swna,a,m = adaptmi(swn;W_m=scale_weighting(cs,params[:W_m_σ]),
                     shape_y = x -> max(0,x),progressbar=progressbar,
                     adapt_params...)

  weight_max = maximum(swna)
  if weight_max > reasonable_response_maximum
    throw(ResponseOverflow(weight_max))
  end

  cs .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs).*im)
  if intermediate_results
    cs,a,m
  else
    cs
  end
end
