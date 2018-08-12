using AuditoryModel
using AxisArrays
const reasonable_response_maximum = 100

function bound(x,min,max)
  y = 1/(1+exp(-8((x - min)/(max - min) - 0.5)))
  miny = 1/(1+exp(-8(-0.5)))
  maxy = 1/(1+exp(-8(0.5)))
  min + (max-min)*(clamp(y,miny,maxy) - miny)/(maxy - miny)
end

function strengths(tracks;kwds...)
  strengths = map_components(tracks;kwds...) do components
    component_means(components)
  end

  AxisArray(hcat(strengths...),axes(strengths,Axis{:time}))
end

function ratio_to_lengths(ratio,threshold=2,min_length=1s)
  percept_lengths(AxisArray(ratios .< threshold,
                            axes(ratios,Axis{:time})),min_length)
end

function component_ratio(tracks;min_length=1s,
                         intermediate_results=false,kwds...)
  map_components(tracks;kwds...) do components
    strengths = sort(component_means(components),rev=true)
    strengths[1] / sum(strengths[2:end])
  end
end

function estimate_bandwidth(sp;threshold=0.25,window=500ms,step=250ms)
  map_windowing(sp,length=window,step=step) do window
    means = mean(window,1)
    peak = maximum(means)
    over = find(means .> threshold*peak)
    maximum(over) - minimum(over) + 1
  end
end

function component_bandwidth_ratio(cs,tracks;min_length=1s,
                                   threshold=0.25,
                                   progressbar=false,window=500ms,
                                   step=250ms)
  crmask = mask(cs,tracks,window=window,step=step)
  spmask = audiospect(crmask,progressbar=progressbar)
  sp = audiospect(cs,progressbar=progressbar)

  fullband = estimate_bandwidth(sp,window=window,step=step,
                               threshold=threshold)
  maskband = estimate_bandwidth(spmask,window=window,step=step,
                               threshold=threshold)

  AxisArray(maskband ./ fullband,axes(fullband,Axis{:time}))
end

function meanabs(A,n)
  result = reducedim((x,y) -> x+abs(y),A,n,zero(real(eltype(A))))
  result ./= size(A,n)
  squeeze(result,n)
end

function findweights(condition,x)
  if condition == :scales
    AxisArray(meanabs(x,axisdim(x,Axis{:freq})) .*
              ustrip.(uconvert.(cycoct,AuditoryModel.scales(x)))',
              axes(x,Axis{:time}), axes(x,Axis{:scale}))
  elseif condition ∈ [:freqs,:track]
    x
  end
end

function remove_key_prefix!(prefix,dict)
  for key in keys(dict)
    if startswith(string(key),prefix)
      dict[Symbol(string(key)[length(prefix)+1:end])] = dict[key]
    end
  end
end

apply_bistable(x,args...;kwds...) = apply_bistable!(deepcopy(x),args...;kwds...)
function apply_bistable!(x,condition,params,settings;
                         interactive=false,
                         intermediate_results=interactive,
                         progressbar=interactive)
  if !(params[:condition] == condition ||
       ((params[:condition] == :scales_track) && (condition ∈ [:scales,:track])))
    return (x,)
  end

  if params[:condition] == :scales_track
    if condition == :scales
      remove_key_prefix!("s_",params)
    elseif condition == :track
      remove_key_prefix!("t_",params)
    end
  end

  noise_params = Dict(
    :τ_σ => params[:τ_σ], :c_σ => params[:c_σ]
  )
  adapt_params = Dict(
    :c_m => params[:c_m], :τ_m => params[:τ_m],
    :c_a => params[:c_a], :τ_a => params[:τ_a],
    :c_x => params[:c_x], :τ_x => params[:τ_x]
  )

  weights = findweights(condition,x)
  weights .= bound.(weights,
                    settings[string(condition)]["bistable"]["input_bound"]...)
  input_weights = copy(weights)

  wn = drift(weights;noise_params...,progressbar=progressbar)
  wna,a,m = adaptmi(wn;W_m=weighting(x,condition,params[:W_m_σ],params[:W_m_c]),
                    shape_y = x -> max(0,x),progressbar=progressbar,
                    adapt_params...)

  freq = settings[string(condition)]["bistable"]["lowpass"]
  order = settings[string(condition)]["bistable"]["lowpass_order"]
  low = digitalfilter(Lowpass(freq;fs=ustrip(1/Δt(wna))),Butterworth(order))
  wna_low = filt!(similar(wna),low,wna)

  if eltype(x) <: Complex
    x .= sqrt.(abs.(x) .* max.(0.0,wna_low)) .* exp.(angle.(x).*im)
  else
    x .*= wna_low
  end

  if intermediate_results
    x,input_weights,wna_low,a,m
  else
    (x,)
  end

end
