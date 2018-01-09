using DataFrames
include("units.jl")
include("tempc.jl")
include("stim.jl")
include("adaptmi.jl")
include("cortmi.jl")
setup_sound(sample_rate=8kHz)

R"library(ggplot2)"
quartz() = R"quartz()"
dir = "../../plots/run_2017_01_09"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect)

x = @> ab(120ms,120ms,1,10,500Hz,6) attenuate(10)

sp = spect(x);
cr = cort(sp);

# GOAL: with MI only, make sure one of the scales wins out


params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=0,τ_a=1s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));


cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales2(cort,mean(abs.(m),[2,4]))
plot_scales2(cort,mean(abs.(cr),[2,4]))
plot_scales2(cort,mean(abs.(cra),[2,4]))

tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[3s])
rplot(tempc,C[3s])
rplot(tempc,[C,Ca],"adapt/mi?" => ["no","yes"])

### okay, now that that's happening, can adaptation give
### us some oscilations???

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales2(cort,mean(abs.(m),[2,4]))

plot_scales2(cort,mean(abs.(cr),[2,4]))
plot_scales2(cort,mean(abs.(cra),[2,4]))

tempc = TCAnalysis(cort,1,1s,method=:pca)
# C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[4s])
rplot(tempc,Ca[2s])

rplot(tempc,[C,Ca],"adapt/mi?" => ["no","yes"])

# NOTE: ratios are clearly off at this point...

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

plot_scales2(cort,mean(abs.(cra),[2,4]))

tempc = TCAnalysis(cort,1,1s,method=:pca)

################################################################################
# how are we going to weight the smaller scales?
# through what justification? just hack it?

# somethink like pink noise...
y = rand(192,62).*0.1.+0.9
crn = cort(y)

weights = mean(abs.(crn),[1,2,4])
weights ./= sum(weights)

# what happens if we weight by this?

# TODO...


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
