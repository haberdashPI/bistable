using Revise
using DataFrames
include("units.jl"); Revise.track("units.jl")
include("stim.jl"); Revise.track("stim.jl")
include("tempc.jl"); Revise.track("tempc.jl")
include("adaptmi.jl"); Revise.track("adaptmi.jl")
include("cortmi.jl"); Revise.track("cortmi.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2017_01_18"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1,10))
tempc = TCAnalysis(cort,1,window=1s,method=:pca,frame_len=10ms)

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
save_plot($(joinpath(dir,"mi_only.pdf")),p,base_aspect_ratio=1.3)
"""

C = tempc(cr);
Ca = tempc(cra);

p = rplot(tempc,Ca[4s])
R"""
p = $p + ggtitle("Mutual-inhibition principle component (at 4 seconds)")
save_plot($(joinpath(dir,"mi_only_pc.pdf")),p,base_aspect_ratio=1.1,ncol=2)
"""

masked, = mask(tempc,Ca[4s],cr,phase=-π/2);
sp = inv(cort,masked,usematlab=false)
p = rplot(spect,sp)

R"""
p = $p + ggtitle("Effect of mutual-inhibition mask (from 4 seconds)")
save_plot($(joinpath(dir,"mi_spect.pdf")),p,base_aspect_ratio=1.3)
"""

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
save_plot($(joinpath(dir,"mi_adapt.pdf")),p,base_aspect_ratio=1.3)
"""

tempc = TCAnalysis(cort,1,window=1s,method=:pca,frame_len=10ms)
Ca = tempc(cra)

## trying out the min filtering mask
t1 = 4.8s
rplot(tempc,Ca[t1])
masked,phase = mask(tempc,Ca[t1],cr[1:400,:,:,:],phase=min_filtering);
sp = inv(cort,masked)
p1 = rplot(spect,sp)

t2 = 10s
rplot(tempc,Ca[t2])
masked, = mask(tempc,Ca[t2],cr[1:400,:,:,:],phase=min_filtering);
sp = inv(cort,masked)
p2 = rplot(spect,sp)

# within this range of scales this apporach can always seperate the sounds
# I could try to expand the range, or I could use a different masking approach...

R"""
p1 = $p1 + ggtitle($("Effect of MI & adaptation mask at $t1"))
p2 = $p2 + ggtitle($("Effect of MI & adaptation mask at $t2"))
p = plot_grid(p1,p2,ncol=2)
# save_plot($(joinpath(dir,"mi_adapt_spect.pdf")),p,base_aspect_ratio=1.3,
          # ncol=2)
"""

p1 = rplot(tempc,Ca[t1])
p2 = rplot(tempc,Ca[t2])


R"""
p1 = $p1 + ggtitle($("MI & adaptation mask at $t1"))
p2 = $p2 + ggtitle($("MI & adaptation mask at $t2"))
p = plot_grid(p1,p2,nrow=2)
save_plot($(joinpath(dir,"mi_adapt_pc.pdf")),p,base_aspect_ratio=1.1,
          nrow=2,ncol=2)
"""

########################################
# let's try using the alternative mask


## trying out the min filtering mask
t1 = 4.8s
rplot(tempc,Ca[t1])
masked = mask2(tempc,Ca[t1],cr[1:400,:,:,:]);
sp = inv(cort,masked)
p1 = rplot(spect,sp)

xl_a = @> ab(120ms,120ms,1,50,500Hz,6,:without_b) normpower amplify(-20)
xl_b = @> ab(120ms,120ms,1,50,500Hz,6,:without_a) normpower amplify(-20)
sp_a = spect(xl_a)[1:400,:]
sp_b = spect(xl_b)[1:400,:]

20 * log10(sqrt(mean(sp .* sp_a)) / sqrt(mean(sp .* sp_b)))

20 * log10(sqrt(mean(sp .* target)) /
         sqrt(mean(abs.(sp.^2 .- sp.*target))))
au_sp = inv(spect,sp,iterations=25)

t2 = 10s
rplot(tempc,Ca[t2])
masked = mask2(tempc,Ca[t2],cr[1:400,:,:,:]);
sp = inv(cort,masked)
p2 = rplot(spect,sp)

# let's just look at the overal seperation

spi = mean_spect(tempc,Ca,cr)

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
