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
dir = "../../plots/run_2017_01_17"
isdir(dir) || mkdir(dir)

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

rplot(tempc,C[3s])

f,phases = fusion_ratio(tempc,C,cr);
snr_a, = object_SNR(tempc,C,cr,spect(x_a));
snr_b, = object_SNR(tempc,C,cr,spect(x_b));

df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])

R"""
library(cowplot)
p = ggplot($df,aes(x=time,y,color=stim)) + geom_line() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  coord_cartesian(ylim=c(-5,5)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using maximum energy phase mask.")

save_plot($(joinpath(dir,"max_phase_SNR.pdf")),p,base_aspect_ratio=1.3)
"""

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


snr_a, = object_SNR(tempc,C,cr,spect(x_a),phase=max_filtering);
snr_b, = object_SNR(tempc,C,cr,spect(x_b),phase=max_filtering);

df = DataFrame(y = [snr_a;snr_b],
               time = repeat(ustrip(times(tempc,cr)),outer=2),
               stim = [fill("A",length(snr_a)); fill("B",length(snr_b))])


R"""
ggplot($df,aes(x=time,y,color=stim)) + geom_line() +
  scale_color_brewer(palette='Set1',name='Stimulus') +
  coord_cartesian(ylim=c(-10,10)) +
  xlab('Time (s)') + ylab('SNR (dB)') +
  ggtitle("SNR when using maximum unmasked energy.")
"""
max_e = max.(snr_a,snr_b)

# thinking....
