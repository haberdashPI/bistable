using DataFrames
include("units.jl"); Revise.track("units.jl")
include("stim.jl"); Revise.track("stim.jl")
include("tempc.jl"); Revise.track("tempc.jl")

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=10,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect,rates=sort([-2.^(-1:0.5:5); 2.^(-1:0.5:5)]))


x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-2)
x = x[0s .. 1s]

cr1 = run_unwrapped(cort,spect(x));
cr2 = run_unwrapped2(cort,spect(x));

quartz(); rplot(cort,real.(cr1[:,:,:,3,:]),scales=[0.5,1.5,4],rates=[-8,-4,-2,2,4,8])
quartz(); rplot(cort,real.(cr2[:,:,:,3,:]),scales=[0.5,1.5,4],rates=[-8,-4,-2,2,4,8])

quartz(); rplot(cort,real.(cr1[:,:,:,8,:]),scales=[0.5,1.5,4],rates=[-8,-4,-2,2,4,8])
quartz(); rplot(cort,real.(cr2[:,:,:,8,:]),scales=[0.5,1.5,4],rates=[-8,-4,-2,2,4,8])

cr = cort(spect(x))
rplot(cort,cr[:,:,:,:],scales=[0.5,1.5,4],rates=[-8,-4,-2,2,4,8])
