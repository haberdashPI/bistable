push!(LOAD_PATH,"packages")
using AuditoryModel

using RCall
quartz() = R"quartz()"
using AxisArrays
using Sounds
include("util/stim.jl")
x = @> ab(120ms,120ms,1,6,500Hz,6) normpower amplify(-20)
y = audiospect(x)
cr = cortical(y,rates=default_rates,scales=default_scales)


# rplot(cr,rates=[-8Hz,-4Hz,-2Hz,2Hz,4Hz,8Hz],scales=[0.5cycoct, 1cycoct, 2cycoct])

cr = cortical(y,rates=default_rates)
rplot(cr,rates=[-8Hz,-4Hz,-2Hz,2Hz,4Hz,8Hz])
yinv = audiospect(cr)
rplot(yinv)


cs = cortical(y,scales=default_scales)
rplot(cs,scales=[0.5cycoct, 1cycoct, 2cycoct])
yinv = audiospect(cs)
rplot(yinv)

# crs1 = cortical(cr,scales=default_scales)
cs = cortical(y,scales=default_scales)
# rplot(cs,scales=[0.5cycoct, 1cycoct, 2cycoct])
crs2 = cortical(cs,rates=default_rates)

# rplot(crs1,rates=[-8Hz,-4Hz,-2Hz,2Hz,4Hz,8Hz],scales=[0.5cycoct, 1cycoct, 2cycoct])
# quartz(); rplot(crs2,rates=[-8Hz,-4Hz,-2Hz,2Hz,4Hz,8Hz],scales=[0.5cycoct, 1cycoct, 2cycoct])

yinv = audiospect(cr)
rplot(yinv)

yinv = audiospect(crs2)
rplot(yinv)

xinv = Sound(y)
rplot(audiospect(xinv))

xinv = Sound(yinv,target_error=0.025,max_iterations=100)

# test named tuples
cparams = Dict(:rates=>default_rates,:scales=>default_scales)
cr = cortical(y;cparams...)
