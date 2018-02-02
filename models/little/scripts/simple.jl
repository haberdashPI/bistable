push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
include("stim.jl")

using DSP
using Unitful: ustrip
using RCall
using DataFrames
delta_t = 10ms

R"""
library(ggplot2)
"""

dir = "../../run_2012_12_13"

function ggplot_qual(y,t;title="",file=nothing)
  ii = CartesianRange(size(y))
  df = DataFrame(y=y[:],x=Float64.(ustrip.(t[map(x -> x[1],ii)[:]])),
                 unit=map(x -> x[2],ii)[:])

R"""
  y = ggplot($df,aes(x,y,color=factor(unit),group=unit)) + geom_line() +
    theme_classic(base_size=14) +
    scale_color_brewer(palette='Set1',name='Unit') +
    ggtitle($title) + xlab('Time (s)') + ylab('response')

  if($(file != nothing)){
    ggsave(paste($dir,$file),width=6,height=5)
  }else{
    y
  }

"""
end

# running basic examples without excittion
t = 0s:delta_t:50s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

y,a,m = adaptmi(2spikes,AdaptMI(τ_y=10ms,
                               c_m=1,τ_m=30ms,
                               c_a=2,τ_a=1s,shape_y=sig,Δt = delta_t))
ggplot_qual(sig.(y[1:1000,:]),t)


y,a,m = adaptmi(2spikes,AdaptMI(τ_y=10ms,
                               c_m=3,τ_m=30ms,
                               c_a=2,τ_a=1s,shape_y=sig,Δt = delta_t))
ggplot_qual(sig.(y[1:1000,:]),t)
ggplot_qual(sig.(y),t)


# introduce a more irregular input

spikey = zeros(Float64,(length(t),4))
step = floor(Int,200ms / delta_t)
spikey[1:step:end,1] = 1
spikey[3:step:end,2] = 1
spikey[5:step:end,3] = 1
spikey[7:step:end,4] = 1
ftype = digitalfilter(Lowpass(20,fs=floor(Int,1s/delta_t)),Butterworth(3))
spikey[:,1] = filt(ftype,spikey[:,1])
spikey[:,2] = filt(ftype,spikey[:,2])
spikey[:,3] = filt(ftype,spikey[:,3])
spikey[:,4] = filt(ftype,spikey[:,4])

ggplot_qual(spikey[1:300,:],t[1:300])

y,a,m = adaptmi(2spikey,AdaptMI(τ_y=10ms,
                               c_m=3,τ_m=30ms,
                               c_a=2,τ_a=1s,shape_y=sig,Δt = delta_t))
ggplot_qual(y[1:300,:],t[1:300])


y,a,e,m = adaptmi(2spikey,AdaptMI(τ_y=10ms,
                                c_m=4,τ_m=30ms,
                                c_a=20,τ_a=1s,shape_y=sig,Δt = delta_t,
                                c_e=3,τ_e=300ms))
ggplot_qual(y[1:1000,:],t[1:1000])

ggplot_qual(e[1:300,:],t[1:300])


y,a,e,m = adaptmi(2spikey,AdaptMI(τ_y=10ms,
                                c_m=4,τ_m=30ms,
                                c_a=20,τ_a=1s,shape_y=sig,Δt = delta_t,
                                c_e=3.6,τ_e=300ms))
ggplot_qual(y[1:1000,:],t[1:1000])

ggplot_qual(e[1:300,:],t[1:300])

# takeaway - the excitation can have a small, but noticable effect
# on the duration of the fluctuations


y,a,e,m = adaptmi(2spikey,AdaptMI(τ_y=10ms,
                                c_m=20,τ_m=30ms,
                                c_a=21,τ_a=1s,shape_y=sig,Δt = delta_t,
                                c_e=1.5,τ_e=70ms))
ggplot_qual(y[1:1000,:],t[1:1000])
ggplot_qual(e[1:300,:],t[1:300])

ggplot_qual(a[1:1000,:],t[1:1000])
ggplot_qual(m[1:1000,:],t[1:1000])

# however, it's not entirely clear whether this approach will
# allow for ongoing oscilations, it seems that when it is strong
# enough, one spike overpowers the others

# I may need to try a different excitation term...

function ggplot_heat(y,t,title="",file=nothing)
  ii = CartesianRange(size(y))
  df = DataFrame(response=y[:],x=Float64.(ustrip.(t[map(x -> x[1],ii)[:]])),
                 unit=map(x -> x[2],ii)[:])
R"""
  y = ggplot($df,aes(x=x,y=unit,fill=response)) + geom_raster() +
    theme_classic(base_size=14) +
    scale_fill_gradient2(low='black',mid='red',high='yellow',
      name='Response',midpoint=0.5) +
    ggtitle($title) + xlab('Time (s)') + ylab('Unit')

  if($(file != nothing)){
    ggsave(paste($dir,$file),width=6,height=5)
  }else{
    y
  }

"""
end
 

## now... let's ramp up the numbers....
t = 0s:delta_t:50s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*20ms .<= t,i] = 1
end

y,a,m = adaptmi(2spikes,AdaptMI(τ_y=25ms,
                               c_m=0.5,τ_m=30ms,
                               c_a=0.6,τ_a=1s,shape_y=sig))

ggplot_heat(y,t)
ggplot_qual(y,t)


y,a,m = adaptmi(2spikes,AdaptMI(τ_y=25ms,
                               c_m=2,τ_m=30ms,
                               c_a=0.6,τ_a=1s,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(2spikes,AdaptMI(τ_y=25ms,
                               c_m=10,τ_m=50ms,
                               c_a=5,τ_a=3s,shape_y=sig))
ggplot_qual(sig.(y),t)

y,a,m = adaptmi(2noise(spikes,100ms,0.3),
                AdaptMI(τ_y=25ms,
                        c_m=10,τ_m=50ms,
                        c_a=5,τ_a=3s,shape_y=sig))
ggplot_qual(sig.(y),t)
ggplot_qual(y,t)


y,a,m = adaptmi(3noise(spikes,100ms,0.3),
                AdaptMI(τ_y=25ms,
                        c_m=2,τ_m=50ms,
                        c_a=1,τ_a=3s,shape_y=sig))
ggplot_qual(sig.(y),t)
ggplot_qual(y,t)
ggplot_heat(sig.(y),t)

k = 5^2
w = [1-exp(-(ii[1] - ii[2])^2/k) for ii in CartesianRange((20,20))]
w ./= sum(w,1)
W_m(x) = w*x

y,a,m = adaptmi(3noise(spikes,100ms,0.3),
                AdaptMI(τ_y=25ms,
                        c_m=2,τ_m=50ms,W_m=W_m,
                        c_a=1,τ_a=3s,shape_y=sig))
ggplot_heat(sig.(y),t)


y,a,m = adaptmi(3spikes,
                AdaptMI(τ_y=25ms,
                        c_m=2,τ_m=50ms,W_m=W_m,
                        c_a=1,τ_a=3s,shape_y=sig))
ggplot_heat(sig.(y),t)
