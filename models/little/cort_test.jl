include("units.jl")
include("tempc.jl")
include("stim.jl")

setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect,rates=[2,8,32,64],scales=[0.125,1,4,16])
tempc = TCAnalysis(cort=cort,sparsity=0.95,frame_len=3,prior=8.0,τ=100ms)

# x = playable(sound("../../test.wav"))[0s .. 1s,:left]
# y = cort(spect(x));
# plot_cort(cort,y)

x = aba(60ms,60ms,20,1kHz,6);
y = cort(spect(x));
yc = tempc(y);

trange = floor(Int,1.5s / Δt(spect)):floor(Int,2.5s / Δt(spect))
rplot(spect,spect(x)[trange,:])

R"quartz()"

trange = floor(Int,1.5s / Δt(tempc)):floor(Int,2.5s / Δt(tempc))
rplot(tempc,yc[trange,:,:,:])

