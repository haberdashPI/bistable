using Pkg
Pkg.activate(joinpath(@__DIR__,".."))

using Unitful
using AxisArrays
using TOML
using RCall
using VegaLite
using DataFramesMeta
using LinearAlgebra
using Statistics
using StatsBase
using Dates
using ProgressMeter
using Colors

const srcdir = @__DIR__
include(joinpath(srcdir,"count_lengths.jl"))
include(joinpath(srcdir,"measures.jl"))
include(joinpath(srcdir,"plotting.jl"))

R"options(rcalljl_options=list(width=800,height=400))"
