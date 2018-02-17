push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_15"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram(len=25,min_freq=250Hz,max_freq=1kHz)
cort = CorticalModel(spect,scales=2.0.^linspace(-0.5,2,6)) # 3 ?
cohere = CoherenceModel(cort,3,window=500ms,method=:nmf,delta=250ms,
                        minwindow=500ms,maxiter=200,tol=1e-3)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6);
cr = cort(sp);
C = cohere(cr);

#=
using JLD
W = C.W
H = C.H
@save "/Users/davidlittle/Data/bistable_cache.jld" W H
=#

#=
using JLD
@load "/Users/davidlittle/Data/bistable_cache.jld"
C = AuditoryCoherence.NMFSeries(W,H,0.025s,0.25s)
=#

rplot(cohere,C)

cohere = CoherenceModel(cort,3,window=500ms,method=:nmf,delta=250ms,
                        minwindow=500ms,maxiter=10_000,tol=1e-4)

x = @>(ab(120ms,120ms,1,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);
C = cohere(cr);

rplot(cohere,C)

spC = mean_spect(cohere,C,cr)
rplot(spect,spC)

spC = mean_spect(cohere,C,cr,component=2)
rplot(spect,spC)

################################
# now with synchronous tones

cohere = CoherenceModel(cort,3,window=500ms,method=:nmf,delta=250ms,
                        minwindow=500ms,maxiter=2_000,tol=1e-3)

x = @>(ab(120ms,120ms,0,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);
C = cohere(cr);

rplot(cohere,C)

spC = mean_spect(cohere,C,cr)
rplot(spect,spC)

spC = mean_spect(cohere,C,cr,component=2)
rplot(spect,spC)
