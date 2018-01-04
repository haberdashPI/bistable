using DataFrames
include("units.jl")
include("tempc.jl")
include("stim.jl")
include("adaptmi.jl")
setup_sound(sample_rate=8kHz)

dir = "../../plots/run_2017_01_02"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect)

x = @> ab(120ms,120ms,1,10,500Hz,6) attenuate(10)

sp = spect(x);
cr = cort(sp);

function scale_weighting(cort)
  s = log.(scales(cort))
  W = @. 1 - exp(-(s - s')^2 / 0.5)
  W ./= sum(W,1)

  function helper(x)
    m = W*vec(mean(x,[1,3]))
    ones(x) .* reshape(m,1,:,1)
  end
end

function plot_scales(cort,m)
  sm = m[:,1,:,1]
  iis = collect(CartesianRange(size(sm)))
  df = DataFrame(level = vec(sm),
                 time = vec(map(x -> ustrip(times(cort,m)[x[1]]),iis)),
                 scale = vec(map(x -> ustrip(scales(cort)[x[2]]),iis)))

R"""
  library(RColorBrewer)
  pal = brewer.pal(5,'Reds')[2:5]
  ggplot($df,aes(x=time,y=level,group=scale,
                 color=factor(scale),linetype=factor(scale))) +
    geom_line(color='black',linetype='solid',size=1.2) + geom_line() +
    scale_color_manual(values = rep(pal,each=3)[1:11]) +
    scale_linetype_manual(values = rep(c("dotdash","longdash","solid"),4)[1:11])
"""
end

params = AdaptMI(c_m=10,τ_m=200ms,W_m=scale_weighting(cort),
                 c_a=3,τ_a=1s,shape_y = x -> 5max(0,x));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

rplot(cort,cra,scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
rplot(cort,a,scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
# rplot(cort,m,scales=[0.5,2,4,8],rates=[-32,-8,-2,2,8,32])
plot_scales(cort,m)

quartz(); plot_scales(cort,abs.(sum(cr,[2,4])))
quartz(); plot_scales(cort,abs.(sum(cra,[2,4])))

# what do the scales look like without mutual inhibition
params = AdaptMI(c_m=0,τ_m=200ms,W_m=scale_weighting(cort),
                 c_a=3,τ_a=1s,shape_y = x -> max(0,x));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

quartz(); plot_scales(cort,mean(abs.(cra),[2,4]))
quartz(); plot_scales(cort,mean(abs.(cr),[2,4]))


# what do the scales look like without adaptation
params = AdaptMI(c_m=5,τ_m=200ms,W_m=scale_weighting(cort),
                 c_a=0,τ_a=1s,shape_y = x -> max(0,x));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

quartz(); plot_scales(cort,mean(abs.(cr),[2,4]))
quartz(); plot_scales(cort,mean(abs.(cra),[2,4]))


params = AdaptMI(c_m=20,τ_m=200ms,W_m=scale_weighting(cort),
                 c_a=0,τ_a=1s,shape_y = x -> max(0,x));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

quartz(); plot_scales(cort,mean(abs.(cr),[2,4]))
quartz(); plot_scales(cort,mean(abs.(cra),[2,4]))


params = AdaptMI(c_m=100,τ_m=200ms,W_m=scale_weighting(cort),
                 c_a=6,τ_a=1s,shape_y = x -> max(0,x));

cra = similar(cr);
cra,a,m = (adaptmi(cra,params) do cr_t,t,dt_cr
  cr[t,:,:,:]
end);

quartz(); plot_scales(cort,mean(abs.(cr),[2,4]))
quartz(); plot_scales(cort,mean(abs.(cra),[2,4]))

tempc = TCAnalysis(cort,1,1s,method=:pca)
C = tempc(cr)
Ca = tempc(cra)

rplot(tempc,[C,Ca],"adapt/mi?" => ["no","yes"])

# none of the changes lead to fusing of the sound...
# in fact in general, a point we can glean from this is
# that in TC this particular signal is always separated...

# I can see two directions from here:

# 1. it's possible that emphasizing the farther end of the empahsized region
# would lead to a single object parsing.  We would need more competition between
# closely neighboring scales for that to work

# okay, that seems to work in cort_test.jl

# 2. it's possible we need to empahsize the larger scales
# in which case we need to rethink how to introduce MI

# another observation, there are different rates of loss, I think
# this is because there is more dispertion of energy at the larger scales
