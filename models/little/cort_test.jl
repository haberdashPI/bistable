# include("units.jl")
include("tempc.jl")
include("stim.jl")

setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect,rates=[2,8,32,64],scales=[0.125,1,4,16])
tempc = TCAnalysis(cort=cort,sparsity=0.8,frame_len=1,prior=8.0,τ=500ms)

# x = playable(sound("../../test.wav"))[0s .. 1s,:left]
# y = cort(spect(x));
# plot_cort(cort,y)

x = aba(60ms,60ms,10,1kHz,6);
y = cort(spect(x));
yc = tempc(y);

rplot(spect,spect(x))
R"quartz()"
rplot(tempc,yc)

better_tempc = TCAnalysis(cort=cort,sparsity=0.1,frame_len=1,prior=8.0,τ=5s)
byc = better_tempc(spect(x))

rplot(spect,spect(x))
R"quartz()"
rplot(better_tempc,byc)


y = spect(x)
a,Y = eigs(y'*y,nev=25)



# TODO: figure out what parts of the analysis
# are screwing up the output
#
# leaky covariance (ccipca can fix)
# sparsity (ccipca makes irrelevant)
# cortical representation
#   - note: probably not, since tempc(spect(x)) doesn't work
