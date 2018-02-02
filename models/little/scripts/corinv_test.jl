push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect,rates=[-8,-2,-1,1,2,8],scales=[0.5,2,8],bandonly=false)
tempc = TCAnalysis(cort,1,window=2s,method=:pca,frame_len=100ms)

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)
sp = spect(x);
cr = cort(sp);
C = tempc(cr);

m = mask(tempc,C[2.15s],cr[1:80,:,:,:])

sp1,debug1 = inv(cort,m,usematlab=true)
sp2,debug2 = inv(cort,m,usematlab=false)

# debug1_ = mat"gen_cort(1,256,40,[1 1]);"
# mat"HR = conj([$debug1_; zeros(256, 1)]);"
# mat"HR = [HR(1); conj(flipud(HR(2:end)))];"
# mat"HR(256+1) = abs(HR(256+2));"

# debug1_ = mat"HR + 0;"

# debug2_ = rate_filter(-1,256,1000 / spect.len,:low,true)
