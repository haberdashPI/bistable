using DataFrames
include("units.jl")
include("tempc.jl")
include("stim.jl")
include("adaptmi.jl")
include("cortmi.jl")
setup_sound(sample_rate=8kHz)

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
# dir = "../../plots/run_2017_01_18"
# isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect)
tempc = TCAnalysis(cort,1,window=2s,method=:pca,frame_len=100ms)

x_a = @> ab(120ms,120ms,1,10,500Hz,6,:without_a) attenuate(20)
x_b = @> ab(120ms,120ms,1,10,500Hz,6,:without_b) attenuate(20)
x = @> ab(120ms,120ms,1,10,500Hz,6) attenuate(20)

sp = spect(x);
cr = cort(sp);
C = tempc(cr);

# f,phases = fusion_ratio(tempc,C,cr);
snr_a = object_SNR2(tempc,C,cr,spect(x_a));
snr_b = object_SNR2(tempc,C,cr,spect(x_b));

df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])

R"""
library(cowplot)
p = ggplot($df,aes(x=time,y,color=stim)) + geom_line() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  # coord_cartesian(ylim=c(-5,5)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using maximum energy phase mask.")

# save_plot($(joinpath(dir,"max_phase_SNR.pdf")),p,base_aspect_ratio=1.3)
"""

sep = ab_match(tempc,C,cr,x_a,x_b)

masked = mask2(tempc,C[21],cr);
rplot(spect,inv(cort,masked))

sp_masked = mean_spect(tempc,C,cr)

rplot(spect,real.(sp_masked))

aud_spm = inv(spect,real.(sp_masked))

# TODO: examine why the SNR is so wacky
# - is the measure inaccurate? what do the resulting
#   spectrograms look like?
# - try using selective scale adaptmi setup,
# does this get less random?

snr_a, = object_SNR(tempc,C,cr,spect(x_a),phase=min_energy);
snr_b, = object_SNR(tempc,C,cr,spect(x_b),phase=min_energy);

df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])

R"""
p = ggplot($df,aes(x=time,y,color=stim)) + geom_line() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  coord_cartesian(ylim=c(-5,5)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using minimum energy phase mask.")
save_plot($(joinpath(dir,"min_phase_SNR.pdf")),p,base_aspect_ratio=1.3)
"""


snr_a, = object_SNR(tempc,C,cr,spect(x_a),phase=π);
snr_b, = object_SNR(tempc,C,cr,spect(x_b),phase=π);

df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])

R"""
p = ggplot($df,aes(x=time,y,color=stim)) + geom_line() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  coord_cartesian(ylim=c(-10,10)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using phase 0 mask.")
save_plot($(joinpath(dir,"fixed_phase_SNR.pdf")),p,base_aspect_ratio=1.3)
"""

# THOUGHT: would it be faster to calculate
# a series of real filters and compute
# principle components from that?
snr_a, = object_SNR(tempc,C,cr,spect(x_a),phase=max_filtering);
snr_b, = object_SNR(tempc,C,cr,spect(x_b),phase=max_filtering);

df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])


R"""
p = ggplot($df,aes(x=time,y,color=stim)) + geom_line() + geom_point() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  coord_cartesian(ylim=c(-10,10)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using maximum unmasked energy.")

save_plot($(joinpath(dir,"max_unmasked_energy_SNR.pdf")),p,base_aspect_ratio=1.3)
"""


# THOUGHT: would it be faster to calculate
# a series of real filters and compute
# principle components from that?
snr_a, = object_SNR(tempc,C,cr,spect(x_a),phase=min_filtering);
snr_b, = object_SNR(tempc,C,cr,spect(x_b),phase=min_filtering);
alert()
df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])


R"""
p = ggplot($df,aes(x=time,y,color=stim)) + geom_line() + geom_point() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  coord_cartesian(ylim=c(-10,10)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using minimum unmasked energy.")

save_plot($(joinpath(dir,"min_unmasked_energy_SNR.pdf")),p,base_aspect_ratio=1.3)
"""

## let's try the real-filter approach
tempc = TCAnalysis(cort,4,window=1s,method=:real_pca,frame_len=500ms)
C = tempc(cr);

rplot(tempc,C[2])
rplot(tempc,C[3])
rplot(tempc,C[4])
rplot(tempc,C[5])
rplot(tempc,C[6])

# can we succefully use this as a mask
# nSOR = normalized scene to object ratio
function nSOR(object,scene)
  object = object ./ maximum(abs.(object))
  scene = scene ./ maximum(abs.(scene))
  -10log10(mean(object.^2) / mean(scene.^2))
end

m = mask(tempc,C[4],cr)
masked = inv(cort,m);
rplot(spect,masked)
nSOR(masked,sp)

m = mask(tempc,C[4],cr,component=2)
masked = inv(cort,m)
rplot(spect,masked)
nSOR(masked,sp)

rplot(spect,sp)

# okay, so, it looks like we get both interpretations, at about
# equal strength. what if we...

# - make df larger?
sp12 = @> ab(120ms,120ms,1,10,500Hz,12) attenuate(20) spect
C12 = tempc(sp12)
rplot(tempc,C12[4])

m = mask(tempc,C12[4],sp12,component=1);
masked = inv(cort,m)
rplot(spect,masked)
nSOR(masked,sp12)

m = mask(tempc,C12[4],sp12,component=3);
masked = inv(cort,m)
rplot(spect,masked)
nSOR(masked,sp12)

# - make df smaller
sp2 = @> ab(120ms,120ms,1,10,500Hz,1) attenuate(20) spect
C2 = tempc(sp2)
rplot(tempc,C2[5])

m = mask(tempc,C2[5],sp2,component=1);
masked = inv(cort,m)
rplot(spect,masked)
nSOR(masked,sp2)

m = mask(tempc,C2[4],sp2,component=2);
masked = inv(cort,m)
rplot(spect,masked)
nSOR(masked,sp2)


# - alter the emphasis on a give scale?


# - alternative implementaiton: is there a way to generate these figures
#   using the complex component? (plot different phases of the compoent)
#   that might make this faster to implement and easier (?) to explain
tempc = TCAnalysis(cort,4,window=2s,method=:pca,frame_len=500ms)
C = tempc(cr);
rplot(tempc,C[4])

pc = reshape(eigvecs(C[4])[:,1],length(scales(cort)),:)
rplot(spect,real.(pc*exp(-π/4*im)))
rplot(spect,real.(pc*exp((π + π/4)*im)))

sum(abs.(real.(pc*exp((π + π/4)*im))))

phases = linspace(-π,π)
mags = [sum(abs.(real.(pc*exp(ϕ*im)))) for ϕ in phases]

rplot(spect,real.(pc*exp(phases[15]*im)))
rplot(spect,real.(pc*exp(phases[38]*im)))

# short answer... I don't think so
