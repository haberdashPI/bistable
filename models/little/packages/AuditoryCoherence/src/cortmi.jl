export scale_weighting, freq_weighting, weighting

function weighting(x,condition,σ,c)
  condition == :scales ? scale_weighting(x,σ,c) :
  condition == :freqs ? freq_weighting(x,σ,c) :
  condition == :track ? track_weighting(x,σ,c) :
  error("No condition named $condition.")
end

function scale_weight(cort,σ,c)
  s = log.(ustrip.(scales(cort)))
  W = @. c*(1 - exp(-(s - s')^2 / (σ*log(2))^2))

  W
end

function scale_weighting(cort,σ,c)
  W = scale_weight(cort,σ,c)
  function helper(x::AbstractArray{T,2}) where T
    reshape(W*vec(sum(x,2)),:,1)
  end

  function helper(x::AbstractArray{T,3}) where T
    reshape(W*vec(sum(x,(1,3))),1,:,1)
  end

  function helper(x::AbstractArray{T,1}) where T
    W*x
  end

  helper
end

function track_weighting(tracks,σ,c)
  f = log.(ustrip.(axisvalues(axes(tracks,Axis{:prior}))[1]))
  W = @. c*(1 - exp(-(f - f')^2 / (σ*log(2))^2))

  x -> W*x
end

function freq_weighting(spect,σ,c)
  f = log.(ustrip.(freqs(spect)))
  W = @. c*(1 - exp(-(f - f')^2 / (σ*log(2))^2))

  x -> W*x
end
