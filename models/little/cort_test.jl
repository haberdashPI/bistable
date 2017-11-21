include("units.jl")
include("tempc.jl")
include("stim.jl")

setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5")
cort = CorticalModel(spect,rates=[2,8,32,64],scales=[0.125,1,4,16])
tempc = TCAnalysis(cort=cort,sparsity=0.96)

# x = playable(sound("../../test.wav"))[0s .. 1s,:left]
# y = cort(spect(x));
# plot_cort(cort,y)

x = aba(60ms,60ms,4,1kHz,6)
y = cort(spect(x))
yc = tempc(y);

plot_cort(cort,y)
R"""
ggsave('../../plots/run_2017_11_21/cort_6st.pdf')
"""

plot_cort(cort,yc)
R"""
ggsave('../../plots/run_2017_11_21/temp_6st.pdf')
"""

x = aba(60ms,4,1kHz,12);
y = cort(spect(x));
yc = tempc(y);

