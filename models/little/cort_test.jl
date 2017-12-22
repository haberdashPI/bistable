include("units.jl")
include("tempc.jl")
include("stim.jl")
setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",
                            len=25)
sts = [0.1,1,6,12]

cort = CorticalModel(spect,scales=2.^(-1:0.5:4),bandonly=false)
dir = "../../plots/run_2017_12_18"

p = Array{Any}(3)

# reduce the time to compute by limiting the range of frequencies
# employed
f_ab(st) = @>(ab(120ms,120ms,1,10,500Hz,st),attenuate(10))
ab_f = [f_ab(st) for st in sts]
crs = [cort(spect(ab_f[i]),usematlab=false) for i in eachindex(ab_f)];

tempc = TCAnalysis(cort,1,1s,method=:pca)
Cs = [tempc(crs[i]) for i in eachindex(ab_f)];

signal = [fusion_signal(tempc,Cs[i],crs[i])
          for i in eachindex(ab_f)]
# p[1] = plot_resps(signal,sts,"delta f (st)")
p[2] = plot_resps(signal,sts,"delta f (st)")

masked = mask(tempc,Cs[4][3s],crs[4],π);

rplot(cort,masked,rates=[2],scales=[2])
rplot(cort,crs[4],rates=[2],scales=[2])
sp = inv(cort,masked,usematlab=true)
rplot(spect,sp)

masked = mask(tempc,Cs[4][3s],crs[4],0);

rplot(cort,masked,rates=[2],scales=[2])
rplot(cort,crs[4],rates=[2],scales=[2])
sp = inv(cort,masked,usematlab=true)
rplot(spect,sp)


## testing inversion
cr #=,Y,HS,z1,z=# = cort(spect(ab_f[4]),usematlab=false);
rplot(cort,cr,rates=[4],scales=[2])
sp = inv(cort,cr,usematlab=true)
rplot(spect,sp)

rplot(cort,cr,rates=[-8,-4,-2,2,4,8],scales=[2,8,16])

## testing inversion: I'm lost here, why isn't this working without matlab?
cr2 #=, Y2, HS2, z12, z2=# = cort(spect(ab_f[4]),usematlab=true);
rplot(cort,cr2,rates=[16],scales=[2])
sp2 = inv(cort,cr2,usematlab=true)
rplot(spect,sp2)

rplot(cort,cr2,rates=[-8,-4,-2,2,4,8],scales=[2,8,16])


# how does weighting the scales go?

center = 32
r_weights = exp.((log.(rates(cort)) .- center).^2 ./ )
wr_cr = cr .* 
