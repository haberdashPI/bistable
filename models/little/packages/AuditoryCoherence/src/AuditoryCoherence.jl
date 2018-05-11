module AuditoryCoherence
using AuditoryModel
using DataFrames
using Requires

import AuditoryModel: Δt, Δf, times, freqs, scales, rates

export adaptmi, drift, scale_weighting, ncomponents, nunits, CoherenceModel,
    fusion_ratio, object_SNR, mask, scene_object_ratio,
    object_SNR2, ab_match, mean_spect, mean_spect2, AdaptMI

using ProgressMeter
next!(x::Progress) = ProgressMeter.next!(x)
next!(x::Void) = nothing

include("adaptmi.jl")
include("cortmi.jl")
include("cohere.jl")
include("nmf.jl")
include("simple_tracking.jl")
include("prior_tracking.jl")

# after we stop revising this package use can use the conditionl
# dependency (it interferes with 'Revise' automatic updates)
@require RCall include(joinpath(@__DIR__,"rplots.jl"))
#  include(joinpath(@__DIR__,"rplots.jl"))
# @require Gadfly include("gplots.jl")

end
