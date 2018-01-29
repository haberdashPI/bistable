using DataFrames
include("../units.jl")
include("../stim.jl")
include("../audio_spect.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)

y = spect(x)
x_inv = inv(spect,y,usematlab=false)

