module AuditoryCoherence
using AuditoryModel
using DataFrames

import AuditoryModel: rplot, Δt, Δf, times, freqs 

export adaptmi, drift, scale_weighting, ncomponents, CoherenceModel,
    fusion_ratio, object_SNR, mask, scene_object_ratio,
    object_SNR2, ab_match, mean_spect, rplot
    
include("online_pca.jl")
include("adaptmi.jl")
include("cortmi.jl")
include("tempc.jl")

end
