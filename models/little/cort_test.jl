include("units.jl")
include("tempc.jl")
include("stim.jl")
setup_sound(sample_rate=8kHz)

# TODO: possible shorten the frame length?
spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)

f_ab_async(st) = @>(ab(60ms,60ms,1,20,500Hz,st),attenuate(10))
ab_async = [f_ab_async(st) for st in [0.1,1,6,12]]

a_times = 2*2*60.*[0; 1:2:18] .* ms .+ 30ms
b_times = a_times .+ 2*60ms .+ 30ms
ap_times = 2*2*60.*[0; 2:2:19] .* ms .+ 30ms

a_indices = max.(1,round(Int,a_times / Δt(spect)))
b_indices = max.(1,round(Int,b_times / Δt(spect)))
ap_indices = max.(1,round(Int,ap_times / Δt(spect)))

dist = ABDist(a_indices,b_indices,ap_indices)

cort = CorticalModel(spect)
tempc = TCAnalysis(cort,10,250ms)

crs = [cort(ab_async[i]) for i in eachindex(ab_async)];
Cs = [tempc(crs[i]) for i in eachindex(ab_async)];

function plot_resps(x,times)
  df = DataFrame(resp = vcat((real.(x[i]) for i in 1:length(x))...),
                 st = repeat([0.1,1,6,12],inner=length(x[1])),
                 time = repeat(times,outer=length(x)))

R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=resp,group=factor(st),color=factor(st))) +
    geom_line() +
    scale_color_brewer(palette='Set1')
"""
end

dir = "../../plots/run_2017_12_12"

signals1 = [fusion_signal(tempc,Cs[i],crs[i]) for i in eachindex(ab_async)]
plot_resps(signals1,times(cort,crs[1]))
R"ggsave(paste($dir,'/pca10_ratio.pdf',sep=''))"

signals2 = [fusion_signal(tempc,Cs[i],crs[i],dist) for i in eachindex(ab_async)]
plot_resps(signals2,ustrip(a_times))
R"ggsave(paste($dir,'/pca10_abdiff.pdf',sep=''))"

tempc_dir = TCAnalysis(cort,10,250ms,method=:direct)
Cds = [tempc_dir(crs[i]) for i in eachindex(ab_async)];
signals3 = [fusion_signal(tempc_dir,Cds[i],crs[i],dist)
            for i in eachindex(ab_async)]
plot_resps(signals3,ustrip(a_times))
R"ggsave(paste($dir,'/Cx_abdiff.pdf',sep=''))"

tempc_s = TCAnalysis(cort,1,250ms)
Cs_s = [tempc_s(crs[i]) for i in eachindex(ab_async)];
signals4 = [fusion_signal(tempc_s,Cs_s[i],crs[i]) for i in eachindex(ab_async)]
signals5 = [fusion_signal(tempc_s,Cs_s[i],crs[i],dist) for i in eachindex(ab_async)]

plot_resps(signals4,times(cort,crs[1]))
R"ggsave(paste($dir,'/online_pca_ratio.pdf',sep=''))"

plot_resps(signals5,ustrip(a_times))
R"ggsave(paste($dir,'/online_pca_abdiff.pdf',sep=''))"

# TODO: find a template for a and b, then use Cs_s, above and generate a
# continuous measure of separation, perhaps using the distnace from a to
# b. Issue - this will not be super robust to changing the task or setup.
