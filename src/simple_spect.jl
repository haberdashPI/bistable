push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_06"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram(len=25)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1,10),bandonly=true)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:real_pca,delta=50ms,
                        n_phases=12)
Cr = cohere(cr);

p1 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6,:without_b)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=:real_pca,delta=50ms,
                        n_phases=12)
Cr = cohere(cr);

p3 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,0,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=:real_pca,delta=50ms,
                        n_phases=12)
Cr = cohere(cr);

p5 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p6 = rplot(spect,sp)

R"""
p = plot_grid($p2 + ggtitle("Alternating AB"),
              $p4 + ggtitle("A alone"),
              $p6 + ggtitle("Synchronized AB"),
              $p1,$p3,$p5,
              nrow=2,ncol=3,rel_heights=c(0.25,0.75))
save_plot($(joinpath(dir,"1_real_valued_components.pdf")),p,
  base_aspect_ratio=3,nrow=5,ncol=3)
"""
########################################
# complex valued components

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,delta=50ms)
Cr = cohere(cr);

p1 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6,:without_b)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=:pca,delta=50ms)
Cr = cohere(cr);

p3 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,0,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=:pca,delta=50ms)
Cr = cohere(cr);

p5 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p6 = rplot(spect,sp)

R"""
p = plot_grid($p2 + ggtitle("Alternating AB"),
              $p4 + ggtitle("A alone"),
              $p6 + ggtitle("Synchronized AB"),
              $p1,$p3,$p5,
              nrow=2,ncol=3,rel_heights=c(0.25,0.75))
save_plot($(joinpath(dir,"2_complex_valued_components.pdf")),p,
  base_aspect_ratio=3,nrow=5,ncol=3)
"""
########################################
# complex valued with simplified scales (real and fake data)

cort = CorticalModel(spect,scales = [1,4],bandonly=true) # 3 ?
cohere = CoherenceModel(cort,6,window=750ms,method=:pca,delta=50ms)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6)
cr = cort(sp);
Cr = cohere(cr);

p1 = rplot(cohere,Cr);
p2 = rplot(spect,sp)

cort = CorticalModel(spect,scales = [1,4],bandonly=true) # 3 ?
cohere = CoherenceModel(cort,6,window=750ms,method=:pca,delta=50ms)

x = ab(120ms,120ms,1,6,500Hz,6) |> normpower |> amplify(-10dB)
sp = spect(x)
cr = cort(sp);
Cr = cohere(cr);

p1 = rplot(cohere,Cr)
p2 = rplot(spect,sp)

spC = mean_spect(cohere,Cr,cr)
rplot(spect,spC)

spC = mean_spect(cohere,Cr,cr,component=2)
rplot(spect,spC)

########################################
# complex vlaued with NMF (fake data)

spect = AuditorySpectrogram(len=25,min_freq=250Hz,max_freq=1kHz)
cort = CorticalModel(spect,scales=2.0.^linspace(-0.5,2,6)) # 3 ?
cohere = CoherenceModel(cort,6,window=500ms,method=:nmf,delta=250ms,
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

########################################
# complex valued components with real stimulis

x = ab(120ms,120ms,1,6,500Hz,6) |> normpower |> amplify(-10dB)
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,delta=50ms)
Cr = cohere(cr);

p1 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

x = ab(120ms,120ms,1,6,500Hz,6,:without_b) |> normpower |> amplify(-10dB)
sp = spect(x);
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,delta=50ms)
Cr = cohere(cr);

p3 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

x = ab(120ms,120ms,0,6,500Hz,6) |> normpower |> amplify(-10dB)
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,delta=50ms)
Cr = cohere(cr);

p5 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p6 = rplot(spect,sp)

R"""
p = plot_grid($p2 + ggtitle("Alternating AB"),
              $p4 + ggtitle("A alone"),
              $p6 + ggtitle("Synchronized AB"),
              $p1,$p3,$p5,
              nrow=2,ncol=3,rel_heights=c(0.25,0.75))
save_plot($(joinpath(dir,"3_realstim_complex_valued_components.pdf")),p,
  base_aspect_ratio=3,nrow=5,ncol=3)
"""



########################################
# real valued components with real stimulis

x = ab(120ms,120ms,1,6,500Hz,6) |> normpower |> amplify(-10dB)
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,9,window=750ms,method=:real_pca,delta=50ms,
                        n_phases=12)
Cr = cohere(cr);

p1 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

x = ab(120ms,120ms,1,6,500Hz,6,:without_b) |> normpower |> amplify(-10dB)
sp = spect(x);
cr = cort(sp);

cohere = CoherenceModel(cort,9,window=750ms,method=:real_pca,delta=50ms,
                        n_phases=12)
Cr = cohere(cr);

p3 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

x = ab(120ms,120ms,0,6,500Hz,6) |> normpower |> amplify(-10dB)
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,9,window=750ms,method=:real_pca,delta=50ms,
                        n_phases=12)
Cr = cohere(cr);

p5 = rplot(cohere,Cr,scales=[0.25,0.5,1,2]);
p6 = rplot(spect,sp)

R"""
p = plot_grid($p2 + ggtitle("Alternating AB"),
              $p4 + ggtitle("A alone"),
              $p6 + ggtitle("Synchronized AB"),
              $p1,$p3,$p5,
              nrow=2,ncol=3,rel_heights=c(0.25,0.75))
save_plot($(joinpath(dir,"4_realstim_real_valued_components.pdf")),p,
  base_aspect_ratio=3,nrow=5,ncol=3)
"""

########################################
# complex masking of simplified stimulus

spect = AuditorySpectrogram(len=25)
cort = CorticalModel(spect,
                     scales = 2.0.^linspace(-2,1.5,10),
                     rates = [-2.0.^(0:0.5:5); 2.0.^(0:0.5:5)],
                     bandonly=false)
cohere = CoherenceModel(cort,4,window=750ms,method=:pca,delta=50ms,
                        normalize_phase=true)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6);
cr = cort(sp);
C = cohere(cr)

rplot(cohere,C,scales=[0.25,0.5,1,2])

spC = mean_spect(cohere,C,cr)
p1 = rplot(spect,spC)

spC = mean_spect(cohere,C,cr,component=2)
p2 = rplot(spect,spC)

R"""
p = plot_grid($p1 + ggtitle("Ideal alternating AB, masking with component 1"),
              $p2 + ggtitle("Ideal alternating AB, masking with component 2"),
              nrow=2,ncol=1)
save_plot($(joinpath(dir,"5_alter_ab_complex_masked.pdf")),p,
  base_aspect_ratio=3,nrow=2,ncol=1)
"""

x = ab(120ms,120ms,1,6,500Hz,6) |> normpower |> amplify(-10dB)
sp = spect(x)
cr = cort(sp);
C = cohere(cr)

rplot(cohere,C,scales=[0.25,0.5,1,2])

spC = mean_spect(cohere,C,cr)
p1 = rplot(spect,spC)

spC = mean_spect(cohere,C,cr,component=2)
p2 = rplot(spect,spC)

C.u[:,:,1] .*= exp(π/2*im)
spC = mean_spect(cohere,C,cr)
p1_ = rplot(spect,spC)

R"""
p = plot_grid($p1 + ggtitle("Alternating AB, masking with component 1"),
              $p2 + ggtitle("Alternating AB, masking with component 2"),
              $p1_ + ggtitle("Alternating AB, masking with component 1 (rotated by pi/2)"),
              nrow=3,ncol=1)
save_plot($(joinpath(dir,"6_realstim_alter_ab_complex_masked.pdf")),p,
  base_aspect_ratio=3,nrow=3,ncol=1)
"""
