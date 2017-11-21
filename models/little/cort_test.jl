include("units.jl")
include("tempc.jl")
include("stim.jl")

setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5")
cort = CorticalModel(spect,rates=2.0.^(1:5),scales=2.0.^(-2:3))
tempc = TCAnalysis(cort=cort,sparsity=0.95)

# x = playable(sound("../../test.wav"))[0s .. 1s,:left]
# y = cort(spect(x));
# plot_cort(cort,y)

x = aba(60ms,4,1kHz,6)
y = cort(spect(x))
yc = tempc(y)

plot_cort(cort,y)
plot_cort(cort,yc)


x = aba(60ms,4,1kHz,12);
y = cort(spect(x));
yc = tempc(y);

plot_cort(cort,y)
R"""
ggsave('../../plots/run_2017_11_16/cort_12st.pdf')
"""

plot_cort(cort,yc)
R"""
ggsave('../../plots/run_2017_11_16/tempc_12st.pdf')
"""


x = aba(60ms,4,1kHz,3);
y = cort(spect(x));
yc = tempc(y);

plot_cort(cort,y)
R"""
ggsave('../../plots/run_2017_11_16/cort_3st.pdf')
"""


plot_cort(cort,yc)
R"""
ggsave('../../plots/run_2017_11_16/tempc_3st.pdf')
"""


x = aba(60ms,100ms,4,1kHz,3);
y = cort(spect(x));
plot_cort(cort,y)
R"""
ggsave('../../plots/run_2017_11_16/cort_wgaps_3st.pdf')
"""


x = aba(60ms,100ms,4,1kHz,12);
y = cort(spect(x));
plot_cort(cort,y)
R"""
ggsave('../../plots/run_2017_11_16/cort_wgaps_12st.pdf')
"""

cort = CorticalModel(spect,
                     rates=2.0.^linspace(1,5,4),
                     scales=2.0.^linspace(-2,3,4))
tempc = TCAnalysis(cort=cort,sparsity=0.95,τ=100ms)

x = aba(60ms,75ms,4,1kHz,3);
y = cort(spect(x))
yc = tempc(y)
