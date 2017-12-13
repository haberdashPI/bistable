include("units.jl")
include("tempc.jl")
include("stim.jl")
setup_sound(sample_rate=8kHz)

# TODO: possible shorten the frame length?
spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=10)
sts = [0.1,1,6,12]

# a_times = 2*2*60.*[0; 1:2:18] .* ms .+ 30ms
# b_times = a_times .+ 2*60ms .+ 30ms
# ap_times = 2*2*60.*[0; 2:2:19] .* ms .+ 30ms

# a_indices = max.(1,round.(Int,a_times / Δt(spect)))
# b_indices = max.(1,round.(Int,b_times / Δt(spect)))
# ap_indices = max.(1,round.(Int,ap_times / Δt(spect)))

# dist = ABDist(a_indices,b_indices,ap_indices)

cort = CorticalModel(spect,scales=2.^(-1:0.5:4))
dir = "../../run_2012_12_13"

function plot_resps(x,vars,varname)
  @show vars
  df = DataFrame(resp = vcat((real.(x[i])
                              for i in 1:length(x))...),
                 var = vcat((fill(vars[i],length(x[i]))
                             for i in eachindex(x))...),
                 time = vcat((ustrip(eachindex(xi) * Δt(spect))
                              for xi in x)...))

p = R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=resp,group=factor(var),color=factor(var))) +
    geom_line() +
    scale_color_brewer(palette='Set1',name=$varname) +
    coord_cartesian(ylim=c(1.1,0)) + ylab('lambda^2 / var(x)')
"""
end

tempc = TCAnalysis(cort,1,1s)

p = Array{Any}(3)

f_ab_async(st) = @>(ab(120ms,120ms,1,10,500Hz,st),attenuate(10))
ab_async = [f_ab_async(st) for st in sts]
crs_async = [cort(ab_async[i]) for i in eachindex(ab_async)];
Cs_async = [tempc(crs_async[i]) for i in eachindex(ab_async)];

sig_async = [fusion_signal(tempc,Cs_async[i],crs_async[i])
             for i in eachindex(ab_async)]
p[1] = plot_resps(sig_async,sts,"delta f (st)");

f_ab_async_d(delta) = @>(ab(delta,delta,1,10,500Hz,12),attenuate(10))
deltas = [60ms, 80ms, 100ms, 120ms]
abd_async = [f_ab_async_d(delta) for delta in deltas]
crd_async = [cort(abd_async[i]) for i in eachindex(abd_async)];
Cd_async = [tempc(crd_async[i]) for i in eachindex(abd_async)];
sigd_async = [fusion_signal(tempc,Cd_async[i],crd_async[i])
              for i in eachindex(abd_async)]

p[2] = plot_resps(sigd_async,2ustrip(deltas),"delta t (ms)");

f_ab_sync(st) = @>(ab(120ms,120ms,0,10,500Hz,st),attenuate(10))
ab_sync = [f_ab_sync(st) for st in sts]
crs_sync = [cort(ab_sync[i]) for i in eachindex(ab_sync)];
Cs_sync = [tempc(crs_sync[i]) for i in eachindex(ab_sync)];
sig_sync = [fusion_signal(tempc,Cs_sync[i],crs_sync[i])
            for i in eachindex(ab_sync)]
p[3] = plot_resps(sig_sync,sts,"delta f (st)");

R"""
library(cowplot)

p = plot_grid($(p[1]) + ggtitle("Asynchronous, by delta f (delta t = 240 ms)"),
          $(p[2]) + ggtitle("Asynchronous, by delta t (delta f = 12 st)"),
          $(p[3]) + ggtitle("Synchronous, by delta f (delta t = 240 ms)"),
          align = "h", ncol = 3)
save_plot($(dir*"/eigenratio.pdf"),p,base_aspect_ratio=1.3,ncol=3,nrow=1)
"""
