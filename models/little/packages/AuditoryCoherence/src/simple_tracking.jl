# steps:

# at each step take an MAP approach
# 1. for each scale:
#   a. eliminate components below a threshold
#   b. group remaing observed components with their closest source
#   c. any source not included is marked as inactive
#   d. use the model prior to compute a score
# 2. select the model with the highest score

# NOTE: not quite gonna work, since NMF doens't work well when there's
# noise.

function track_sources(C::NMFSeries,params::SimpleTracker)
end
