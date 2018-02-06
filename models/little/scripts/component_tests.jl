push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_01"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram(len=25,min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1,10))

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)

sp = spect(x);
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=(:real_pca,8),frame_len=50ms)
Cr = cohere(cr);
scale_plot(cohere,Cr,scales=[0.25,0.5,1,2])

cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=50ms)
Cc = cohere(cr);
quartz(); scale_plot(cohere,Cc,scales=[0.25,0.5,1,2])

scaleonly = CorticalModel(spect,scales = scales(cort),rates=[NaN])
quartz(); rplot(scaleonly,scaleonly(sp),scales=[0.25,0.5,1,2])

