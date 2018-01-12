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
dir = "../../plots/run_2017_01_09"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect)
tempc = TCAnalysis(cort,1,window=1s,method=:pca,frame_len=100ms)

x_a = @> ab(120ms,120ms,1,10,500Hz,6,:without_a) attenuate(20)
x_b = @> ab(120ms,120ms,1,10,500Hz,6,:without_b) attenuate(20)
x = mix(x_a,x_b)

sp = spect(x);
cr = cort(sp);
C = tempc(cr);

rplot(tempc,C[3s])

m = mask(tempc,C,cr);
masked = inv(cort,m,usematlab=true);
p = rplot(spect,masked)



scene_SNR(x,a,b) = 20 * abs(log10(sqrt(mean(a .* x)) / sqrt(mean(b .* x))))

# IDEA: find the maimum amplitude element in the component and
# use its phase, apply the mask at that time, and compute SNR
# (also would be worth plotting the resulting spectrogram,
#  changing the component as we go)


# GOAL: with MI only, make sure one of the scales wins out

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
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

masked = mask(tempc,Ca[4s],cr,0);
sp = inv(cort,masked,usematlab=true)
p = rplot(spect,sp)

R"""
p = $p + ggtitle("Effect of mutual-inhibition mask (from 4 seconds)")
save_plot($(joinpath(dir,"mi_spect.pdf")),p,base_aspect_ratio=1.3)
"""

### okay, now that that's happening, can adaptation give
### us some oscilations???

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales2(cort,mean(abs.(cra),[2,4]),name="mean |x|")

########################################
# let's try a longer run of that...

xl = @> ab(120ms,120ms,1,50,500Hz,6) attenuate(10)

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cr = cort(xl);

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra),[2,4]),name="mean |x|")
R"""
p = $p + ggtitle("MI & adaptation of scales")
save_plot($(joinpath(dir,"mi_adapt.pdf")),p,base_aspect_ratio=1.3)
"""

tempc = TCAnalysis(cort,1,1s,method=:pca)
Ca = tempc(cra)

rplot(tempc,Ca[4.2s])
masked = mask(tempc,Ca[4.2s],cr[1:400,:,:,:],0);
sp = inv(cort,masked,usematlab=true)
p1 = rplot(spect,sp)


rplot(tempc,Ca[8s])
masked = mask(tempc,Ca[8s],cr[1:400,:,:,:],π);
sp = inv(cort,masked,usematlab=true)
p2 = rplot(spect,sp)

R"""
p1 = $p1 + ggtitle("Effect of MI & adaptation mask at 4.2s")
p2 = $p2 + ggtitle("Effect of MI & adaptation mask at 8s")
p = plot_grid(p1,p2,ncol=2)
save_plot($(joinpath(dir,"mi_adapt_spect.pdf")),p,base_aspect_ratio=1.3,
          ncol=2)
"""

p1 = rplot(tempc,Ca[4.2s])
p2 = rplot(tempc,Ca[8s])


R"""
p1 = $p1 + ggtitle("MI & adaptation mask at 4.2s")
p2 = $p2 + ggtitle("MI & adaptation mask at 8s")
p = plot_grid(p1,p2,nrow=2)
save_plot($(joinpath(dir,"mi_adapt_pc.pdf")),p,base_aspect_ratio=1.1,
          nrow=2,ncol=2)
"""

################################################################################
# how are we going to weight the smaller scales?
# through what justification? just hack it?

# somethink like pink noise...
y = rand(192,62).*0.1.+0.9

p = rplot(spect,y)
R"""
save_plot($(joinpath(dir,"noisy_spect.pdf")),
          $p + ggtitle("Noisy spectrogram"),
          base_aspect_ratio=1.3)
"""

crn = cort(y)

weights = mean(abs.(crn),[1,2,4])
weights ./= sum(weights)

R"""
p = qplot(x=$(scales(cort)),y=$(1./vec(weights)),geom='line') +
  ggtitle("Scale weights derived from noisy spectrogram") +
  xlab("Scale (cyc/oct)") + ylab("Weight")
save_plot($(joinpath(dir,"scale_weights.pdf")),p,base_aspect_ratio=1.3)
"""

# what happens if we weight by this?

# let's just look at mutual inhibition and if something still wins out,
# (and what is it?)

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=0,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cr = cort(x);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]./weights[1,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra),[2,4]))
R"""
p = $p + ggtitle("Scale weighted MI")
save_plot($(joinpath(dir,"mi_only_scale_weight.pdf")),p,base_aspect_ratio=1.3)
"""

# tempc = TCAnalysis(cort,1,1s,method=:pca)
# C = tempc(cr);
Ca = tempc(cra);

masked = mask(tempc,Ca[4s],cr,-π/2);
sp = inv(cort,masked,usematlab=true)
p = rplot(spect,sp)

R"""
p = $p + ggtitle("Effect of scale weighted MI")
save_plot($(joinpath(dir,"mi_only_scale_weight_spect.pdf")),p,
          base_aspect_ratio=1.3)
"""

########################################

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cr = cort(x);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]./weights[1,:,:,:]
end);

p = plot_scales2(cort,mean(abs.(cra),[2,4]))
R"""
p = $p + ggtitle("Scale weighted Adapt/MI")
save_plot($(joinpath(dir,"mi_adapt_scale_weight.pdf")),p,base_aspect_ratio=1.3)
"""

p = plot_scales2(cort,mean(abs.(cra[100:end,:,:,:]),[2,4]))
R"""
p = $p + ggtitle("Scale weighted Adapt/MI (second half)")
save_plot($(joinpath(dir,"mi_adapt_scale_weight_half2.pdf")),p,base_aspect_ratio=1.3)
"""


tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[3s])
rplot(tempc,C[3s])

# mmm... k, not super convincing, try it over a longer time period

cr = cort(xl);
cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]./weights[1,:,:,:]
end);

plot_scales2(cort,mean(abs.(m),[2,4]))
plot_scales2(cort,mean(abs.(a),[2,4]))
plot_scales2(cort,mean(abs.(cr./weights),[2,4]))
plot_scales2(cort,mean(abs.(cra),[2,4]))
plot_scales2(cort,mean(abs.(cra[400:600,:,:,:]),[2,4]))

tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[3s])
rplot(tempc,C[3s])

# nope... I think what's happening is that whatever property leads to poor
# winning out of a single reponse also prevents the faster rates from winning
# out over the slower rates. This is probably related to the comment in ranking
# about need an extra term to handle the gaps between tones



# xl = @> ab(120ms,120ms,1,100,500Hz,6) attenuate(10)
# spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
#                             min_freq = 250Hz,max_freq=1500Hz)
# cort = CorticalModel(spect)

# params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
#                  c_a=8,τ_a=2s,shape_y = x -> max(0,x),
#                  Δt = Δt(cort));

# cr = cort(xl);
# cr ./= weights;

# cra = similar(cr);
# cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
#   cr[t,:,:,:]
# end);

# plot_scales(cort,mean(abs.(cra[100:200,:,:,:]),[2,4]))
# rplot(spect,squeeze(mean(abs.(cra[100:end,:,:,:]),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(cra[100:200,:,:,:]),[2,4]),(2,4)))


# # things seem to drop too quickly in this case...
# # let's try decreasing constants

# params = AdaptMI(c_m=2,τ_m=250ms,W_m=scale_weighting(cort,0.5),
#                  c_a=0.5,τ_a=2s,shape_y = x -> max(0,x),
#                  c_σ = 1,τ_σ = 100ms,
#                  Δt = Δt(cort));

# cra = similar(cr);
# cra25,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
#   0.25cr[t,:,:,:]
# end);

# cra100,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
#   cr[t,:,:,:]
# end);

# # rplot(spect,squeeze(mean(abs.(m),[2,4]),(2,4)))
# # rplot(spect,squeeze(mean(abs.(a),[2,4]),(2,4)))
# # rplot(spect,squeeze(mean(abs.(cra),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(cra25[100:400,:,:,:]),[2,4]),(2,4)))
# quartz(); rplot(spect,squeeze(mean(abs.(cra100[100:400,:,:,:]),[2,4]),(2,4)))

# rplot(spect,squeeze(mean(abs.(cra[100:end,:,:,:]),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(cra[100:200,:,:,:]),[2,4]),(2,4)))

# # what if we introduce some noise?

# cra = similar(cr);
# crn = drift(cr,params);

# # debugging noise function
# rplot(spect,abs.(crn[:,18,7,:]))

# cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
#   3crn[t,:,:,:]
# end);

# rplot(spect,squeeze(mean(abs.(crn[100:200,:,:,:]),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(cra),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(cra[100:400,:,:,:]),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(cra[100:200,:,:,:]),[2,4]),(2,4)))


# rplot(spect,squeeze(mean(abs.(m[100:400,:,:,:]),[2,4]),(2,4)))
# rplot(spect,squeeze(mean(abs.(a[100:400,:,:,:]),[2,4]),(2,4)))
