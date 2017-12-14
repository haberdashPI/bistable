using Unitful: ustrip, s, ms
# using Plots; plotlyjs()
using RCall
using DataFrames
include("units.jl")
include("adaptmi.jl")
Δt = 10ms

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

t = 0s:Δt:50s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1


y,a,m = adaptmi(spikes,AdaptMI(τ_y=50ms,
                               c_m=0.2,τ_m=200ms,
                               c_a=0.5,τ_a=3s,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(sig.(y),t)
ggplot_qual(a,t)


t = 0s:Δt:200s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

y,a,m = adaptmi(2spikes,AdaptMI(τ_y=50ms,
                               c_m=0.2,τ_m=300ms,
                               c_a=0.5,τ_a=10s,shape_y=sig))
ggplot_qual(sig.(y),t)
ggplot_qual(a,t)
ggplot_qual(m,t)

y,a,m = adaptmi(2spikes,AdaptMI(τ_y=50ms,
                               c_m=0.5,τ_m=300ms,
                               c_a=0.6,τ_a=5s,shape_y=sig))
ggplot_qual(sig.(y),t)



y,a,m = adaptmi(2spikes,AdaptMI(τ_y=25ms,
                               c_m=0.5,τ_m=300ms,
                               c_a=0.6,τ_a=5s,shape_y=sig))
ggplot_qual(sig.(y),t)
ggplot_qual(a,t)
ggplot_qual(m,t)


y,a,m = adaptmi(2spikes,AdaptMI(τ_y=25ms,
                               c_m=0.5,τ_m=150ms,
                               c_a=0.6,τ_a=3s,shape_y=sig))
ggplot_qual(sig.(y),t)
ggplot_qual(a,t)
ggplot_qual(m,t)


y,a,m = adaptmi(2spikes,AdaptMI(τ_y=25ms,
                               c_m=0.5,τ_m=30ms,
                               c_a=0.6,τ_a=1s,shape_y=sig))
ggplot_qual(sig.(y),t)
ggplot_qual(a,t)
ggplot_qual(m,t)


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
t = 0s:Δt:50s
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
