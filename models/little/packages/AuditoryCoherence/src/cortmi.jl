
# the key is getting this weighting matrix right,
# I think I will be best able to do this
# by thinking in terms of the lower dim space first,
# how can this essentially be a competition between different
# scales (based on the sums in each scale, or something like that)
function scale_weighting(cort,σ=1)
  s = log.(scales(cort))
  W = @. 1 - exp(-(s - s')^2 / (σ^2*log(2)))
  W ./= sum(W,2)

  @show W

  function helper(x::AbstractArray{T,3}) where T
    m = W*vec(sum(x,(1,3)))
    x./max.(1e-10,sum(x,(1,3))) .* reshape(m,1,:,1)
    # ones(x) .* reshape(m,1,:,1)
    # TODO: can the above be expressed as a matrix multiplication?
    # I don't think so... I need another abstraction to expression
    # mutual-inhibition in the paper/grant
  end
end
