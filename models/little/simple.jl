using Unitful: ustrip, s, ms
# using Plots; plotlyjs()
using RCall
using DataFrames
include("units.jl")
include("simple_model.jl")
const Δt = 1ms

R"""
library(ggplot2)
"""
dir = "../../plots/meeting_2017_12_07B/"

function ggplot_qual(y,t,title="",file=nothing)
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

t = 0s:Δt:5s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 2
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

# hold on, what's a reasonable time constant for y
# given the eigenspace model? probably about 1s

y,a,m = adaptmi(spikes,AdaptMI(τ_y=1s,c_m=0,c_a=0,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(spikes,AdaptMI(τ_y=50ms,c_m=0,c_a=0.2,τ_a=50s,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)
ggplot_qual(m,t)

y,a,m = adaptmi(spikes,AdaptMI(τ_y=10ms,c_m=5,c_a=0.2,τ_a=5s,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)
ggplot_qual(m,t)


################################################################################
# old

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0,c_a=1,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=1ms,c_m=0,c_a=10,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0,c_a=1,τ_a=5*1.5s,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0,c_a=2,τ_a=5*1.5s,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,c_a=1,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=25ms,c_m=0.2,c_a=1,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)



y,a,m = adaptmi(3spikes,AdaptMI(τ_y=1ms,c_m=0.2,c_a=1,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=1ms,c_m=0.2,c_a=5,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=1ms,c_m=2,c_a=5,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)


y,a,m = adaptmi(3spikes,AdaptMI(τ_y=1ms,c_m=10,c_a=5,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)


y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,c_a=0.1,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,c_a=0.4,shape_y=sig))
ggplot_qual(y,t)


y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,c_a=0.5,shape_y=sig))
ggplot_qual(y,t)


y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,τ_a=4s,c_a=0.55,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,τ_a=4s,c_a=0.8,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,τ_a=4s,c_a=0.8,c_m=0,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.2,τ_a=4s,c_a=1.2,c_m=0,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.0,τ_a=4s,c_a=0.95,c_m=0,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.0,τ_a=4s,c_a=1.2,c_m=0,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.0,τ_a=4s,c_a=0.90,c_m=0,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.0,c_a=0.8,c_m=0,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=5ms,c_m=0.0,c_a=0.8,c_m=1,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(10spikes,AdaptMI(τ_y=50ms,c_m=0.0,c_a=0.8,c_m=1,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)


t = 0s:Δt:5s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[10ms .<= t,2] = 1
spikes[20ms .<= t,3] = 1
spikes[30ms .<= t,4] = 1

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=0,c_a=0.8,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=0,c_a=1.2,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=0,c_a=4,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=0,c_a=20,shape_y=sig))
ggplot_qual(y,t)


y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=.5,c_a=20,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=20,c_a=20,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=20,τ_m=200ms,c_a=20,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=10,τ_m=200ms,c_a=10,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=3,τ_m=200ms,c_a=3,shape_y=sig))
ggplot_qual(y,t)

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=0.5,τ_m=200ms,c_a=0.01,
                                τ_a=5s,shape_y=sig))
ggplot_qual(y,t)
ggplot_qual(a,t)
ggplot_qual(m,t)

t = 0s:Δt:30s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

y,a,m = adaptmi(3spikes,AdaptMI(τ_y=50ms,c_m=0,c_a=1,shape_y=sig))
ggplot_qual(y,t)


t = 0s:Δt:5s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=10,shape=sig)
ggplot_qual(y,t,"4 Unit, MI = 10 (noisy)","simple_4unit_mi_10_noise.pdf")


function ggplot_heat(y,t,title,file)
  ii = CartesianRange(size(y))
  df = DataFrame(response=y[:],x=Float64.(ustrip.(t[map(x -> x[1],ii)[:]])),
                 unit=map(x -> x[2],ii)[:])
R"""
  ggplot($df,aes(x=x,y=unit,fill=response)) + geom_raster() +
    theme_classic(base_size=14) +
    scale_fill_gradient2(low='black',mid='red',high='yellow',
      name='Response',midpoint=0.5) +
    ggtitle($title) + xlab('Time (s)') + ylab('Unit')
  ggsave(paste($dir,$file),width=6,height=5)
"""
end

## now... let's ramp up the numbers....
t = 0s:Δt:2s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig)
ggplot_heat(y,t,"20 Unit, MI = 10","simple_20unit_mi_10.pdf")

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)
ggplot_heat(y,t,"20 Unit, MI = 10 (noise)","simple_20unit_mi_10_noise.pdf")


# what do things look when only some of the inputs have activity?
t = 0s:Δt:2s
spikes = zeros(length(t),20)
for i in 1:10
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)
ggplot_heat(y,t,"20 Unit (half off), MI = 10 (noise)",
            "simple_20unit_halfoff_mi_10_noise.pdf")

## okay, now let's try a case where inhibition is neighbor based
const k = 5^2
W_mi = [1-exp(-(ii[1] - ii[2])^2/k) for ii in CartesianRange((20,20))]

t = 0s:Δt:2s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, MI = 10 (noise), distance-based MI",
            "simple_20unit_dist_mi_10_noise.pdf")



y,a,mi = adaptmi(identity,3spikes,
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, MI = 10, distance-based MI",
            "simple_20unit_dist_mi_10.pdf")


# how does this look over a longer time period?

t = 0s:Δt:20s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,3spikes,
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, MI = 10, distance-based MI (extended)",
            "simple_20unit_dist_mi_10_extended.pdf")



y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, MI = 10 (noise), distance-based MI (extended)",
            "simple_20unit_dist_mi_10_noise_extended.pdf")

# what about when all inputs turn on simultaneously?
spikes = ones(length(t),20)

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, immediate-on MI = 10 (noise), distance-based MI",
            "simple_20unit_immediate_dist_mi_10_noise.pdf")


y,a,mi = adaptmi(identity,3spikes,
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, immediate-on MI = 10, distance-based MI",
            "simple_20unit_immediate_dist_mi_10.pdf")


spikes = zeros(length(t),20)
spikes[:,1:4] = ones(spikes[:,1:4])
spikes[:,15:19] = ones(spikes[:,15:19])
y,a,mi = adaptmi(identity,3spikes,
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, two spots MI = 10, distance-based MI",
            "simple_20unit_two_spot_dist_mi_10.pdf")


spikes = zeros(length(t),20)
spikes[:,1:4] = ones(spikes[:,1:4])
spikes[:,15:19] = ones(spikes[:,15:19])
y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)
ggplot_heat(y,t,"20 Unit, two spots MI = 10 (noise), distance-based MI",
            "simple_20unit_two_spot_dist_mi_10_nosie.pdf")


#=
# what about when we use lateral, rather than global inhibition?
# TODO: below is unfinished, I have no clear conclusion yet
ψ(x,y,k=4) = let t = x-y
    2 / (sqrt(3k) * π^(1/4)) *
        (1 - ((x-y)/k)^2) * exp(-t^2/(2k^2))
end
W_mi = [-3ψ(ii[1],ii[2]) for ii in CartesianRange((20,20))]

t = 0s:Δt:20s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,3spikes,
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)

spikes = ones(length(t),20)

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)


y,a,mi = adaptmi(identity,noise(3spikes,100ms,1.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)
=#

# TODO: look at different subsets of neurons firing
# TODO: find some way to have the lateral inhibition look
# less "saturated"
# TODO: might be worth looking at lateral inhibition setup
# with one or two neurons, since it also has self excitation


# concluding thoughts - I seem to be getting somewhere with this
# I think once I've explored this a little further I need to
# starting trying this out on the full auditory model
# in layer one, to see how this plays out

# TODO: start identifying multiple tracks of work
# I need to be making progress on with this project
