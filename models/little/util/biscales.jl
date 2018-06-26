using AuditoryModel
using AxisArrays
const reasonable_response_maximum = 100

function meanabs(A,n)
  result = reducedim((x,y) -> x+abs(y),A,n,zero(real(eltype(A))))
  result ./= size(A,n)
  squeeze(result,n)
end

bound(x,min,max) = 1/(1+exp(-4((x - min)/(max - min) - 0.5)))

function bistable_scales(cs,params,settings;
                         interactive=false,intermediate_results=interactive,
                         progressbar=interactive)
  if params[:condition] != :scales
    return (cs,)
  end

  noise_params = Dict(
    :τ_σ => params[:τ_σ], :c_σ => params[:c_σ]
  )
  adapt_params = Dict(
    :c_m => params[:c_m], :τ_m => params[:τ_m],
    :c_a => params[:c_a], :τ_a => params[:τ_a],
    :c_x => params[:c_x], :τ_x => params[:τ_x]
  )

  scale_weights = AxisArray(meanabs(cs,axisdim(cs,Axis{:freq})),
                            axes(cs,Axis{:time}), axes(cs,Axis{:scale}))
  scale_weights .= bound.(scale_weights,0.0,
                          settings["scales"]["bistable"]["input_bound"])

  swn = drift(scale_weights;noise_params...,progressbar=progressbar)
  swna,a,m = adaptmi(swn;W_m=scale_weighting(cs,params[:W_m_σ],params[:W_m_c]),
                     shape_y = x -> max(0,x),progressbar=progressbar,
                     adapt_params...)

  freq = settings["scales"]["bistable"]["lowpass"]
  order = settings["scales"]["bistable"]["lowpass_order"]
  low = digitalfilter(Lowpass(freq;fs=ustrip(1/Δt(swna))),Butterworth(order))
  swna_low = filt!(similar(swna),low,swna)

  cs .= sqrt.(abs.(cs) .* max.(0.0,swna_low)) .* exp.(angle.(cs).*im)

  if intermediate_results
    cs,swna_low,a,m
  else
    (cs,)
  end

end
