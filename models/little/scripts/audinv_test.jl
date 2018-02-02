push!(LOAD_PATH,"packages")
using DataFrames
using AuditoryModel
include("stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"

spect = AuditorySpectrogram(len=25)

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)

y = spect(x)
x_inv = inv(spect,y,usematlab=false)

