using StatsBase: rle
using Interpolations

ecdf(x) = ecdf!(copy(x))
function ecdf!(x)
  sort!(x)
  ux, lens = rle(x)
  y = (cumsum(lens).-0.5) ./ length(x)
  fn = interpolate((unique!(x),),y,Gridded(Linear()))
  extrapolate(fn,Line())
end

"""
    KLestimate(p,q)

Given two 1d sets of samples `p` and `q`, estimate the KL divergence
between `P` and `Q`, where `p` was drawn from distribution `P`
and `q` from `Q`.
"""
function KLestimate(p,q;eps=1e-8)
  P,Q = ecdf(p), ecdf(q)
  mean(@.(log( (P(p)-P(p-eps))/
               (Q(p)-Q(p-eps)) )))
end

#=

Problem: we want to have measure that is relative to human error,
but we don't really have the human error for the histogram plot.
If we assume it is log-normal, we can get the mean human error...
is this a reasonable measure? Can I justify this to others.
I really want to just measure the KL div directly,
but how do I place this on some reasonable scale?
I have no measures of individual listeners...

we could get the KL-divergence for both measures,
and simply add the divergence together...
I worry this will be a very different scale though...

since we don't have individual data, we could bootstrap the KL-divergence
of the human data

=#

