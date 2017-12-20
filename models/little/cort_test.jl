# include("units.jl")
include("tempc.jl")
include("stim.jl")
setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",
                            len=25,min_freq = 400.0Hz,max_freq = 2kHz)
sts = [0.1,1,6,12]

cort = CorticalModel(spect,scales=2.^(-1:0.5:4))
dir = "../../plots/run_2017_12_18"

p = Array{Any}(3)

# reduce the time to compute by limiting the range of frequencies
# employed
f_ab(st) = @>(ab(120ms,120ms,1,10,500Hz,st),attenuate(10))
ab_f = [f_ab(st) for st in sts]
crs = [cort(spect(ab_f[i]),false) for i in eachindex(ab_f)];

tempc = TCAnalysis(cort,1,1s,method=:pca)
Cs = [tempc(crs[i]) for i in eachindex(ab_f)];

signal = [fusion_signal(tempc,Cs[i],crs[i])
          for i in eachindex(ab_f)]
p[1] = plot_resps(signal,sts,"delta f (st)")
