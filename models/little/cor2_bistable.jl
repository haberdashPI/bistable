using DataFrames
include("units.jl")
include("tempc.jl")
include("stim.jl")
include("adaptmi.jl")
include("cortmi.jl")
setup_sound(sample_rate=8kHz)

R"library(ggplot2)"
quartz() = R"quartz()"
dir = "../../plots/run_2017_01_02"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect)

x = @> ab(120ms,120ms,1,10,500Hz,6) attenuate(10)

sp = spect(x);
cr = cort(sp);

function plot_scales(cort,m,range=nothing)
  sm = m[:,1,:,1]
  iis = collect(CartesianRange(size(sm)))
  df = DataFrame(level = vec(sm),
                 time = vec(map(x -> ustrip(times(cort,m)[x[1]]),iis)),
                 scale = vec(map(x -> ustrip(scales(cort)[x[2]]),iis)))

  @show unique(df[:scale])
R"""
  library(RColorBrewer)
  pal = brewer.pal(5,'Reds')[2:5]
  p = ggplot($df,aes(x=time,y=level,group=scale,
                 color=factor(round(scale,1)),linetype=factor(round(scale,1)))) +
    geom_line(color='black',linetype='solid',size=1.2) + geom_line() +
    scale_color_manual(values = rep(pal,each=3)[1:11],name="Scale") +
    scale_linetype_manual(values =
      rep(c("dotdash","longdash","solid"),4)[1:11],name="Scale")
"""
  if range != nothing
R"""
      p = p + coord_cartesian(ylim=c($(first(range)),$(last(range))))
"""
  end

  R"p"
end

# none of the changes lead to fusing of the sound...
# in fact in general, a point we can glean from this is
# that in TC this particular signal is always separated...

# I can see two directions from here:

# 1. it's possible that emphasizing the farther end of the empahsized region
# would lead to a single object parsing.  We would need more competition between
# closely neighboring scales for that to work

params = AdaptMI(c_m=5,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=3,τ_a=1s,shape_y = x -> max(0,x));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales(cort,mean(abs.(m),[2,4]))

p = [
  plot_scales(cort,mean(abs.(cr),[2,4]))
  plot_scales(cort,mean(abs.(cra),[2,4]))
];

rplot(spect,squeeze(maximum(abs.(cra),(2,4)),(2,4)))

R"""
library(cowplot)
p = plot_grid($(p[1]) + ggtitle("Without Adapt/MI"),
              $(p[2]) + ggtitle("With Adapt/MI"),align="h",ncol=2)
# save_plot($(joinpath(dir,"adaptmi.png")),p,
#   base_aspect_ratio=1.4,ncol=2)
"""

rplot(cort,m,scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])

tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,[C,Ca],"adapt/mi?" => ["no","yes"])

rplot(tempc,Ca[3s])
rplot(tempc,C[3s])

# 2. it's possible we need to empahsize the larger scales
# in which case we need to rethink how to introduce MI

# another observation, there are different rates of loss, I think
# this is because there is more dispertion of energy at the larger scales


# GOAL: with MI only, make sure one of the scales wins out


params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=0,τ_a=1s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));


cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales(cort,mean(abs.(m),[2,4]))

plot_scales(cort,mean(abs.(cr),[2,4]))
plot_scales(cort,mean(abs.(cra),[2,4]))
 
rplot(cort,m[150:end,:,:,:],scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
rplot(cort,cr[150:end,:,:,:],scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
rplot(cort,cra[150:end,:,:,:],scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])

tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[3s])
rplot(tempc,C[3s])

rplot(tempc,[C,Ca],"adapt/mi?" => ["no","yes"])


### okay, now that that's happening, can adaptation give
### us some osscilations???

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));


cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales(cort,mean(abs.(m),[2,4]))

plot_scales(cort,mean(abs.(cr),[2,4]))
plot_scales(cort,mean(abs.(cra),[2,4]))
 
rplot(cort,m[150:end,:,:,:],scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
rplot(cort,cr[150:end,:,:,:],scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
rplot(cort,cra[150:end,:,:,:],scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])

tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[3s])
rplot(tempc,Ca[4.5s])
rplot(tempc,C[3s])

rplot(tempc,[C,Ca],"adapt/mi?" => ["no","yes"])

# NOTE: ratios are clearly off at this point...


# let's try a longer run of that...


xl = @> ab(120ms,120ms,1,100,500Hz,6) attenuate(10)
spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect)

params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort,0.1),
                 c_a=8,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(cort));

cr = cort(xl);

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

plot_scales(cort,mean(abs.(cra),[2,4]))

tempc = TCAnalysis(cort,1,1s,method=:pca)
# C = tempc(cr);
Ca = tempc(cra);

rplot(tempc,Ca[4.5s])
rplot(tempc,Ca[7s])

p = Array{Any}(2)

m4_5 = mask(tempc,Ca[4.5s],cr,0);
sp4_5 = inv(cort,m4_5,usematlab=true)
p[1] = rplot(spect,sp4_5)

m7 = mask(tempc,Ca[7s],cr,0);
sp7 = inv(cort,m7,usematlab=true)
p[2] = rplot(spect,sp7)

rplot(tempc,[Ca],"adapt/mi?" => ["yes"])
