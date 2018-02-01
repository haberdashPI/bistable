using DataFrames
include("units.jl"); Revise.track("units.jl")
include("stim.jl"); Revise.track("stim.jl")
include("tempc.jl"); Revise.track("tempc.jl")
include("adaptmi.jl"); Revise.track("adaptmi.jl")
include("cortmi.jl"); Revise.track("cortmi.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_01"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1,10))
tempc = TCAnalysis(cort,4,window=750ms,method=(:real_pca,8),frame_len=600ms)

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)

sp = spect(x);
cr = cort(sp);
C = tempc(cr);

# GOAL: with MI only, make sure one of the scales wins out

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=0,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra),[2,4]),name="mean |x|")
R"""
p = $p + ggtitle("Mutual inhibition of scales")
# save_plot($(joinpath(dir,"mi_only.pdf")),p,base_aspect_ratio=1.3)
"""

C = tempc(cr);
Ca = tempc(cra);

p = rplot(tempc,Ca[4s])
R"""
p = $p + ggtitle("Mutual-inhibition principle component (at 4 seconds)")
# save_plot($(joinpath(dir,"mi_only_pc.pdf")),p,base_aspect_ratio=1.1,ncol=2)
"""

sp_C = mean_spect(tempc,C,cr,component=1)
sp_Ca = mean_spect(tempc,Ca,cr,component=1)

# TODO: plot

rplot(spect,sp_C)
rplot(spect,sp_Ca)

### okay, now that that's happening, can adaptation give
### us some oscilations???

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=12,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales2(cort,mean(abs.(cra),[2,4]),name="mean |x|")

########################################
# let's try a longer run of that...

xl = @> ab(120ms,120ms,1,50,500Hz,6) normpower amplify(-10)

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cr = cort(xl);

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra[1:end,:,:,:]),[2,4]),name="mean |x|")
R"""
p = $p + ggtitle("MI & adaptation of scales")
# save_plot($(joinpath(dir,"mi_adapt.pdf")),p,base_aspect_ratio=1.3)
"""

# IDEA: the speed of this needs to slow down!

tempc = TCAnalysis(cort,4,window=750ms,method=(:real_pca,8),frame_len=600ms)
Ca = tempc(cra)

## let's see how this works at the extremes
t1 = 3.8s
t2 = 9.5s
rplot(tempc,Ca[t2])

sp1 = mean_spect(tempc,Ca,cr,component=1)

p1 = plot_scales2(cort,mean(abs.(cra[1:end,:,:,:]),[2,4]),name="mean |x|")
p2 = rplot(tempc,Ca[t1],n=2,showvar=false)
p3 = rplot(tempc,Ca[t2],n=2,showvar=false)
p4 = rplot(spect,sp1)

R"""
p = plot_grid($p2 + ggtitle($("Components at $t1")),
              $p3 + ggtitle($("Components at $t2")),
              $p1 + ggtitle("Mean response for each scale"),
              $p4 + ggtitle("Windowed Masking"),
              nrow=4,ncol=1)

save_plot($(joinpath(dir,"1_bistable_scales.pdf")),p,
  base_aspect_ratio=2,nrow=4,ncol=1)
"""

################################################################################
# OLD STUFF

# IDEA: plot an individual scale slice of a component across time to get
# a sense of how it changes
# scalev(tempc,Ci,scalei) =
#   reshape(eigvecs(Ci)[:,1],length(scales(tempc.cort)),:)[scalei,:]

# scale1_c1 = hcat(scalev.(tempc,Ca,1,8)...)'
# rplot(spect,scale1_c1)

# scale1_c2 = hcat(scalev.(tempc,Ca,2,8)...)'
# quartz(); rplot(spect,scale1_c2)

# next steps:

#=
- look at bistable rates
- look at bistable principal components

- what about role of attention?
- what about role of prediction/wm?

=#

########################################
# let's look at noise

xl = @> ab(120ms,120ms,1,50,500Hz,6) normpower amplify(-10)

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 τ_σ = 100ms,c_σ = 0.045,
                 Δt = Δt(cort));

cr = cort(xl);
crn = drift(cr,params);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  crn[t,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra[100:end,:,:,:]),[2,4]),name="mean |x|")

R"""
p = $p + ggtitle("MI & adaptation & noise (sd = 0.045) over scales")
save_plot($(joinpath(dir,"mi_adapt_noise_045.pdf")),p,base_aspect_ratio=1.6)
"""


params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 τ_σ = 100ms,c_σ = 0.06,
                 Δt = Δt(cort));

crn = drift(cr,params);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  crn[t,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra[100:end,:,:,:]),[2,4]),name="mean |x|")

R"""
p = $p + ggtitle("MI & adaptation & noise (sd = 0.06) over scales")
save_plot($(joinpath(dir,"mi_adapt_noise_06.pdf")),p,base_aspect_ratio=1.6)
"""
