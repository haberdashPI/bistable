push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "plots"
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,2,50Hz,6,:ramp) |> normpower |> amplify(-10dB)

R"""
qplot(x=$((1:length(x))./ustrip(samplerate())),y=$(Array(x)),geom='line') +
  xlab('Time (s)')
ggsave($(joinpath(dir,"timeamp.pdf")))
"""

x = ab(120ms,120ms,1,4,500Hz,6,:ramp) |> normpower |> amplify(-10dB)

y = audiospect(x)
p = rplot(y)
R"""
ggsave($(joinpath(dir,"spect.pdf")))
"""

cs = cortical(y,scales=cycoct.*round.(2.0.^linspace(-2,2,3),1))
rplot(cs,fn=abs)
R"""
ggsave($(joinpath(dir,"scales.pdf")))
"""

csr = cortical(y,scales=cycoct.*round.(2.0.^linspace(-2,2,9),1),
               rates=[Hz.*.-round.(2.0.^linspace(0,5,11),1);
                      Hz.*round.(2.0.^linspace(0,5,11),1)])

rplot(csr[:,[14,18,20],1:3:end,:],fn=abs)
R"""
ggsave($(joinpath(dir,"scale_rates.pdf")))
"""
