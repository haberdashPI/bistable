using AuditoryModel
using AxisArrays
const reasonable_response_maximum = 100

function bound(x,min,max)
  y = 1/(1+exp(-8((x - min)/(max - min) - 0.5)))
  miny = 1/(1+exp(-8(-0.5)))
  maxy = 1/(1+exp(-8(0.5)))
  min + (max-min)*(clamp(y,miny,maxy) - miny)/(maxy - miny)
end

function count_streams(tracks;window=500ms,step=250ms,threshold=2,min_length=1s,
                       progressbar=false,intermediate_results=false)
  Ct1 = tracks[1][1]
  @assert axisdim(Ct1,Axis{:time}) == 1
  windows = windowing(Ct1,length=window,step=step)
  ts = linspace(times(Ct1)[1],times(Ct1)[end]-step,length(windows))
  ratios = fill(0.0,length(windows))

  for (i,ixs) in enumerate(windows)
    best_track = map(tracks) do results
      mean(results[2][ixs])
    end |> indmax
    track_window = tracks[best_track][1][Axis{:time}(ixs)]

    strengths = sort(component_means(track_window),rev=true)
    ratios[i] = strengths[1] / sum(strengths[2:end])
  end

  counts = percept_lengths(AxisArray(ratios .< threshold,Axis{:time}(ts)),
                           min_length)
  if intermediate_results
    counts,AxisArray(ratios,Axis{:time}(ts))
  else
    (counts,)
  end
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

function apply_bistable(x,condition,params,settings;
                        interactive=false,
                        intermediate_results=interactive,
                        progressbar=interactive)
  if params[:condition] != condition
    return (x,)
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
