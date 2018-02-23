push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_20"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram(len=25,min_freq=250Hz,max_freq=1kHz)
cort = CorticalModel(spect,scales=2.0.^linspace(-0.5,2,6)) # 3 ?
cohere = CoherenceModel(cort,3,window=100ms,method=:nmf,delta=50ms,
                        maxiter=200,tol=1e-3)

x = @>(ab(120ms,120ms,1,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);
C = cohere(cr);

p1 = rplot(cohere,C)

spC = mean_spect2(cohere,C,cr);
p2 = rplot(spect,spC)

spC = mean_spect2(cohere,C,cr,component=2);
p3 = rplot(spect,spC)

x = @>(ab(120ms,120ms,0,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);
C = cohere(cr);

p4 = rplot(cohere,C)

spC = mean_spect2(cohere,C,cr)
p5 = rplot(spect,spC)

R"""
p = plot_grid($p1 + ggtitle("Asycnhronous Components"),
              plot_grid($p2 + ggtitle("Masking for Asycnhronous Component 1"),
                        $p3 + ggtitle("Masking for Asycnhronous Component 2"),
                        nrow=1,ncol=2,align='h'),
              $p4 + ggtitle("Synchronized Components"),
              $p5 + ggtitle("Masking for Synchronized Component 1"),
              rel_heights=c(2,1,2,1),nrow=4,ncol=1)
save_plot($(joinpath(dir,"1_nmf_components.pdf")),p,base_aspect_ratio=1.4,nrow=4,ncol=2)
"""
