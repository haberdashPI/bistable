module AuditoryCoherence
using AuditoryModel
using DataFrames
using Requires

import AuditoryModel: Δt, Δf, times, freqs, scales, rates

export adaptmi, drift, scale_weighting, ncomponents, CoherenceModel,
    fusion_ratio, object_SNR, mask, scene_object_ratio,
    object_SNR2, ab_match, mean_spect
    
include("online_pca.jl")
include("adaptmi.jl")
include("cortmi.jl")
include("tempc.jl")

# after we stop revising this package use can use the conditionl
# dependency (it interferes with 'Revise' automatic updates)
# @require RCall include(joinpath(@__DIR__,"rplots.jl"))
include(joinpath(@__DIR__,"rplots.jl"))
# @require Gadfly include("gplots.jl")

end
