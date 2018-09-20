using Pkg
using Unitful
using AxisArrays
Pkg.activate("..")
using TOML
using RCall
const srcdir = @__DIR__
include(joinpath(srcdir,"count_lengths.jl"))

R"options(rcalljl_options=list(width=800,height=400))"
