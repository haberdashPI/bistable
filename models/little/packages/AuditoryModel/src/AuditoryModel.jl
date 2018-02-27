module AuditoryModel
using DataFrames
using Unitful
using Unitful: ms, s, Hz, kHz

export ms, s, Hz, kHz

include("modelresult.jl")
include("audiospect.jl")
include("cortical.jl")

# @require RCall include(joinpath(@__DIR__,"rplots.jl"))
include(joinpath(@__DIR__,"rplots.jl"))

end
