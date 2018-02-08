# push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
# include("stim.jl")
include("../stim.jl")


R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_06"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram(len=25)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1,10),bandonly=true)

function ideal_ab(spect,tone_len,spacing_len,offset_ratio,repeats,freq,delta,options...)
  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a_freqi = indmin(abs.(freqs(spect) .- a_freq))
  b_freqi = indmin(abs.(freqs(spect) .- b_freq))

  a_times = 2(tone_len+spacing_len) .* (0:10)
  b_times = 2(tone_len+spacing_len) .* (0:10) .+
    offset_ratio*(tone_len+spacing_len)
  sp = zeros(ceil(Int,(maximum(b_times) + (tone_len+spacing_len)) / Δt(spect)),
             length(freqs(spect)))

  a_timei = find(mapslices(any,a_times .< times(spect,sp)' .<
                           (a_times .+ (tone_len)),[1]))
  b_timei = find(mapslices(any,b_times .< times(spect,sp)' .<
                           (b_times .+ (tone_len)),[1]))

  if !(:without_a in options) sp[a_timei,a_freqi] = 1 end
  if !(:without_b in options) sp[b_timei,b_freqi] = 1 end

  sp
end

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=(:real_pca,12),frame_len=50ms)
Cr = cohere(cr);

p1 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6,:without_b)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=(:real_pca,12),frame_len=50ms)
Cr = cohere(cr);

p3 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,0,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=(:real_pca,12),frame_len=50ms)
Cr = cohere(cr);

p5 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
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

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,frame_len=50ms)
Cr = cohere(cr);

p1 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6,:without_b)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=50ms)
Cr = cohere(cr);

p3 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

sp = ideal_ab(spect,120ms,120ms,0,6,500Hz,6)
cr = cort(sp);

cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=50ms)
Cr = cohere(cr);

p5 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
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
# complex valued components with real stimulis

x = @>(ab(120ms,120ms,1,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,frame_len=50ms)
Cr = cohere(cr);

p1 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

x = @>(ab(120ms,120ms,1,6,500Hz,6,:without_b),normpower,amplify(-10))
sp = spect(x);
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,frame_len=50ms)
Cr = cohere(cr);

p3 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

x = @>(ab(120ms,120ms,0,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,6,window=750ms,method=:pca,frame_len=50ms)
Cr = cohere(cr);

p5 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
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

x = @>(ab(120ms,120ms,1,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,9,window=750ms,method=(:real_pca,12),frame_len=50ms)
Cr = cohere(cr);

p1 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p2 = rplot(spect,sp)

x = @>(ab(120ms,120ms,1,6,500Hz,6,:without_b),normpower,amplify(-10))
sp = spect(x);
cr = cort(sp);

cohere = CoherenceModel(cort,9,window=750ms,method=(:real_pca,12),frame_len=50ms)
Cr = cohere(cr);

p3 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
p4 = rplot(spect,sp)

x = @>(ab(120ms,120ms,0,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);

cohere = CoherenceModel(cort,9,window=750ms,method=(:real_pca,12),frame_len=50ms)
Cr = cohere(cr);

p5 = scale_plot(cohere,Cr,scales=[0.25,0.5,1,2]);
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
cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=50ms,
                        normalize_phase=true)

sp = ideal_ab(spect,120ms,120ms,1,6,500Hz,6);
cr = cort(sp);
C = cohere(cr)

scale_plot(cohere,C,scales=[0.25,0.5,1,2])

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

x = @>(ab(120ms,120ms,1,6,500Hz,6),normpower,amplify(-10))
sp = spect(x)
cr = cort(sp);
C = cohere(cr)

scale_plot(cohere,C,scales=[0.25,0.5,1,2])

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
