push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_06"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram(len=25,min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1.5,10))
cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=600ms)

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)

sp = spect(x);
cr = cort(sp);
C = cohere(cr);

# GOAL: with MI only, make sure one of the scales wins out

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=0,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

p = collapsed_scale_plot(cort,mean(abs.(cra),[2,4]),name="mean |x|")
R"""
p = $p + ggtitle("Mutual inhibition of scales")
# save_plot($(joinpath(dir,"mi_only.pdf")),p,base_aspect_ratio=1.3)
"""

C = cohere(cr);
Ca = cohere(cra);

p = rplot(cohere,Ca[4s])
R"""
p = $p + ggtitle("Mutual-inhibition principle component (at 4 seconds)")
# save_plot($(joinpath(dir,"mi_only_pc.pdf")),p,base_aspect_ratio=1.1,ncol=2)
"""

sp_C = mean_spect(cohere,C,cr,component=1)
sp_Ca = mean_spect(cohere,Ca,cr,component=1)

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

collapsed_scale_plot(cort,mean(abs.(cra),[2,4]),name="mean |x|")

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

p1 = collapsed_scale_plot(cort,mean(abs.(cra[1:end,:,:,:]),[2,4]),name="mean |x|")

cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=600ms)
Ca = cohere(cra)

p2 = scale_plot(cohere,Ca,scales=[0.25,0.5,1,2])

spC = mean_spect(cohere,Ca,cr)
p3 = rplot(spect,spC)

Ca_ = deepcopy(Ca)
Ca_.u[:,:,1] .*= exp(π/2*im)
spC_ = mean_spect(cohere,Ca_,cr)
p4 = rplot(spect,spC_)

R"""
p = plot_grid($p1 + ggtitle("Mean absolute response for each scale."),
              $p2 + ggtitle("Principle Components by scale"),
              $p3 + ggtitle("Masked response (component 1)."),
              $p4 + ggtitle("Masked response (component 1 rotated by pi/2)"),
              nrow=4,ncol=1)

save_plot($(joinpath(dir,"7_bistable_scales.pdf")),p,
  base_aspect_ratio=2,nrow=4,ncol=1)
"""

################################################################################
# OLD STUFF

# IDEA: plot an individual scale slice of a component across time to get
# a sense of how it changes
# scalev(cohere,Ci,scalei) =
#   reshape(eigvecs(Ci)[:,1],length(scales(cohere.cort)),:)[scalei,:]

# scale1_c1 = hcat(scalev.(cohere,Ca,1,8)...)'
# rplot(spect,scale1_c1)

# scale1_c2 = hcat(scalev.(cohere,Ca,2,8)...)'
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
                 τ_σ = 100ms,c_σ = 0.02,
                 Δt = Δt(cort));

cr = cort(xl);
crn = drift(cr,params);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  crn[t,:,:,:]
end);

p = collapsed_scale_plot(cort,mean(abs.(cra[100:end,:,:,:]),[2,4]),name="mean |x|")

R"""
p = $p + ggtitle("MI & adaptation & noise (sd = $(parasm.c_σ)) over scales")
save_plot($(joinpath(dir,"mi_adapt_noise_045.pdf")),p,base_aspect_ratio=1.6)
"""

Ca = cohere(cra)

p1 = scale_plot(cohere,C,scales=[0.25,0.5,1,2]);

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 τ_σ = 100ms,c_σ = 0.06,
                 Δt = Δt(cort));

crn = drift(cr,params);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  crn[t,:,:,:]
end);

p = collapsed_scale_plot(cort,mean(abs.(cra[100:end,:,:,:]),[2,4]),name="mean |x|")

R"""
p = $p + ggtitle("MI & adaptation & noise (sd = 0.06) over scales")
save_plot($(joinpath(dir,"mi_adapt_noise_06.pdf")),p,base_aspect_ratio=1.6)
"""
