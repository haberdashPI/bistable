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
# cohere = CoherenceModel(cort,4,window=750ms,method=:pca,frame_len=600ms)
cohere = CoherenceModel(cort,3,window=100ms,method=:nmf,delta=50ms,
                        maxiter=200,tol=1e-3)

x = @> ab(120ms,120ms,1,6,500Hz,6) normpower amplify(-20)

sp = spect(x);
cr = cort(sp);
C = cohere(cr);

rplot(cohere,C)

# GOAL: with MI only, make sure one of the scales wins out

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.5),
                 c_a=0,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

p = collapsed_scale_plot(cort,mean(abs.(cra),[2,4]))
R"""
p = $p + ggtitle("Mutual inhibition of scales")
# save_plot($(joinpath(dir,"mi_only.pdf")),p,base_aspect_ratio=1.3)
"""

C = cohere(cr);
Ca = cohere(cra);

rplot(cohere,Ca)

sp_C = mean_spect2(cohere,C,cr,component=1)
sp_Ca = mean_spect2(cohere,Ca,cr,component=1)

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

collapsed_scale_plot(cort,mean(abs.(cra),[2,4]))

########################################
# let's try a longer run of that...

xl = @> ab(120ms,120ms,1,50,500Hz,6) normpower amplify(-10)
cort = CorticalModel(spect,scales=2.0.^linspace(-1.5,2,8))

params = AdaptMI(c_m=20,τ_m=500ms,W_m=scale_weighting(cort,0.5),
                 c_a=8,τ_a=5s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cr = cort(xl);

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

p1 = collapsed_scale_plot(cort,mean(abs.(cra[1:end,:,:,:]),[2,4]))

cohere = CoherenceModel(cort,3,window=100ms,method=:nmf,delta=50ms,
                        maxiter=200,tol=1e-3)
Ca = cohere(cra)

p2 = rplot(cohere,Ca)

spC = mean_spect2(cohere,Ca,cr,component=1)
p3 = rplot(spect,spC)

# spC = mean_spect2(cohere,Ca,cr,component=2)
# p3 = rplot(spect,spC)

# Ca_ = deepcopy(Ca)
# Ca_.u[:,:,1] .*= exp(π/2*im)
# spC_ = mean_spect(cohere,Ca_,cr)
# p4 = rplot(spect,spC_)

R"""
p = plot_grid($p1 + ggtitle("Mean absolute response for each scale."),
              $p2 + ggtitle("Components by scale"),
              $p3 + ggtitle("Masked response (component 1)."),
              nrow=3,ncol=1)

save_plot($(joinpath(dir,"4_bistable_scales.pdf")),p,
  base_aspect_ratio=2,nrow=4,ncol=1)
"""
########################################
# looking at noise

xl = @> ab(120ms,120ms,1,150,500Hz,6) normpower amplify(-10)
params = AdaptMI(c_m=20,τ_m=500ms,W_m=scale_weighting(cort,0.5),
                 c_a=8,τ_a=5s,shape_y = x -> max(0,x),
                 τ_σ = 100ms,c_σ = 0.09,
                 Δt = Δt(cort));

cr = cort(xl);
crn = drift(cr,params);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  crn[t,:,:,:]
end);

p1 = collapsed_scale_plot(cort,mean(abs.(cra[1:end,:,:,:]),[2,4]))

# Thought: the noise as written probably averages out to somethign very
# minimal across all components of a scale, I should probably have
# per scale noise
