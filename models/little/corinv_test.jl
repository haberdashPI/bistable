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

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect,rates=[-8,-2,-1,1,2,8],scales=[0.5,2,8],bandonly=false)
tempc = TCAnalysis(cort,1,window=2s,method=:pca,frame_len=100ms)

x = @> ab(120ms,120ms,1,10,500Hz,6) attenuate(20)
sp = spect(x);
cr = cort(sp);

cr = cort(sp,usematlab=false)

sp1 = inv(cort,cr,usematlab=true)
sp2 = inv(cort,cr,usematlab=false)

debug1_ = mat"gen_cort(1,256,40,[1 1]);"
mat"HR = conj([$debug1_; zeros(256, 1)]);"
mat"HR = [HR(1); conj(flipud(HR(2:end)))];"
mat"HR(256+1) = abs(HR(256+2));"

debug1_ = mat"HR + 0;"

debug2_ = rate_filter(-1,256,1000 / spect.len,:low,true)
