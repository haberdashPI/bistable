using Unitful: ustrip, s, ms
using Plots; plotlyjs()
include("simple_model.jl")
const Δt = 1ms

t = 0s:Δt:2s
spikes = zeros(length(t),4)

spikes[t .% 0.5s .< 0.1s,1] = 1
spikes[t .% 0.5s .< 0.2s,2] = 1
spikes[t .% 0.5s .< 0.3s,3] = 1
spikes[t .% 0.5s .< 0.4s,4] = 1

plot(spikes)

y,a = adaptmi_old(identity,spikes[:,4],300ms,5,0)
y,a = adaptmi_old(identity,spikes[:,4],100ms,5,0)
y,a = adaptmi_old(identity,spikes[:,4],100ms,1,0)

Δt = (1/1000)s
t = 0s:Δt:2s
spikes = zeros(length(t),2)
spikes[t .% 0.5s .< 0.3s,1] = 1
spikes[0.1s .<= (t .% 0.5s) .< 0.4s,2] = 1

y,a = adaptmi_old(identity,spikes,300ms,5,0)
y,a = adaptmi_old(identity,spikes,300ms,5,0.6)


spikes = zeros(length(t),2)
spikes[:,1] = 3
spikes[t .> 10ms,2] = 3

y,a = adaptmi_old(identity,spikes,300ms,5,0.6)
y,a = adaptmi_old(identity,spikes,50ms,5,0.6)
y,a = adaptmi_old(identity,spikes,600ms,5,0.6)

y,a = adaptmi_old(identity,spikes,600ms,5,1.9)

y,a = adaptmi_old(identity,spikes,1.5s,8,1.9)
y,a = adaptmi_old(identity,spikes,1.5s,8,1.9,sig)
y,a = adaptmi_old(identity,spikes,100ms,8,1.9,sig)
y,a = adaptmi_old(identity,spikes,3s,8,1.9,sig)


y,a = adaptmi_old(identity,1.2noise(spikes,100ms,0.1),500ms,8,1.9,sig)
y,a = adaptmi_old(identity,1.2noise(spikes,100ms,0.1),500ms,8,1.2,sig)
y,a = adaptmi_old(identity,1.2noise(spikes,300ms,0.1),2.2s,6,1.9,sig)

y,a = adaptmi_old(identity,[1.2 0.4].*noise(spikes,300ms,0.1),2.2s,6,1.9,sig)

y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,8,1.9,sig)
y,a = adaptmi_old(identity,2noise(spikes,100ms,0.1),1.5s,8,1.9,sig)


y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,8,1.9)

y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),2s,1,0.6)

y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,1,1.9)

y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,0.4,1.9)


# POSSIBLE TODO: final - explore parameter space?? (adjust MI and adaptation?)
# do we make sure that this can follow behavioral dynamics of bi-stability??
# NOTE: show that each element of the model is necessary to get cycles


Δt = s * 1/1000

t = 0s:Δt:2s
spikes = zeros(length(t),3)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1

y,a = adaptmi_old(identity,spikes,1.5s,8,0.2,sig)
y,a = adaptmi_old(identity,spikes,1.5s,8,0.6,sig)
y,a = adaptmi_old(identity,spikes,1.5s,8,1.2,sig)
y,a = adaptmi_old(identity,spikes,1.5s,8,1.6,sig)


t = 0s:Δt:2s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

y,a = adaptmi_old(identity,spikes,1.5s,8,0.2,sig)
y,a = adaptmi_old(identity,spikes,1.5s,8,0.6,sig)
y,a = adaptmi_old(identity,spikes,1.5s,8,1.2,sig)

y,a = adaptmi_old(identity,spikes,0.5s,3,0.6,sig)
y,a = adaptmi_old(identity,spikes,0.5s,3,0.8,sig)
y,a = adaptmi_old(identity,spikes,0.5s,3,0.9,sig)
y,a = adaptmi_old(identity,spikes,0.5s,3,0.95,sig)
y,a = adaptmi_old(identity,spikes,0.5s,3,0.975,sig)
y,a = adaptmi_old(identity,spikes,0.5s,3,1.0,sig)
y,a = adaptmi_old(identity,spikes,0.5s,3,1.2,sig)

# NOTE: that doesn't seem to be working well...
# what if we add noise?

y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,8,0.2,sig)
y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,8,0.6,sig)
y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),1.5s,8,1.2,sig)

y,a = adaptmi_old(identity,noise(spikes,100ms,0.1),0.3s,8,1.2,sig)
y,a = adaptmi_old(identity,noise(4spikes,100ms,0.5),1.5s,8,1.2,sig)
y,a = adaptmi_old(identity,noise(4spikes,100ms,0.1),1.5s,8,0.8,sig)

# NOTE: modeling problem - we still run into funky oscillations, I think there
# is a stiffness that is making it difficult to represent what would happen
# discretely. However, I think these depend on very sensitive conditions that
# would occur with a noisy system.

# what if only two of the inputs are strong?

y,a = adaptmi_old(identity,[1 1 0.2 0.2].*spikes,1.5s,8,1.2,sig)
y,a = adaptmi_old(identity,[1 1 0.2 0.2].*spikes,1.5s,8,1.4,sig)
y,a = adaptmi_old(identity,[1 1 0.2 0.2].*spikes,1.5s,8,1.6,sig)
y,a = adaptmi_old(identity,[1 1 0.2 0.2].*spikes,3s,8,1.6,sig)

# the basic problem is that if all of the inputs are close
# to one another, with a high enough MI
# all of the inputs can be tamped down, and suddenly
# there is no MI, so they all jump up again (together)
# so they all come back on, and so forth.

# FIXED: inserting a delay into mutual inhibition
# seems to help a lot

y,a,mi = adaptmi(identity,spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=0.2,shape=sig)
y,a,mi = adaptmi(identity,spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=1.2,shape=sig)
y,a,mi = adaptmi(identity,spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=2,shape=sig)
y,a,mi = adaptmi(identity,spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=10,shape=sig)
y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=10,shape=sig)

y,a,mi = adaptmi(identity,noise(3spikes,τ_a=100ms,c_a=0.5),τ_mi=1.5s,c_mi=8,
                 shape=50ms,10,sig)

y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=2)
y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=3)

y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=3,shape=sig)

# does this now sustain, absent the noise?
t = 0s:Δt:10s
spikes = zeros(length(t),4)
spikes[000ms .<= t,1] = 1
spikes[100ms .<= t,2] = 1
spikes[200ms .<= t,3] = 1
spikes[300ms .<= t,4] = 1

y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=8,τ_mi=50ms,c_mi=3,shape=sig)

# what does noise do?
y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=8,τ_mi=50ms,
                 c_mi=10,shape=sig)
y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)

# what does this look like when I start all of the inputs at the same time?
t = 0s:Δt:10s
spikes = ones(length(t),4)
y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)
y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=750ms,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)

## now... let's ramp up the numbers....
t = 0s:Δt:2s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)
y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig)

# what do things look when only some of the inputs have activity?
t = 0s:Δt:2s
spikes = zeros(length(t),20)
for i in 1:10
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=10,shape=sig)
y,a,mi = adaptmi(identity,3spikes,τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig)

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),τ_a=1.5s,c_a=3,τ_mi=50ms,
                 c_mi=20,shape=sig)

## okay, now let's try a case where inhibition is neighbor based
const k = 5^2
W_mi = [1-exp(-(ii[1] - ii[2])^2/k) for ii in CartesianRange((20,20))]

t = 0s:Δt:2s
spikes = zeros(length(t),20)
for i in 1:20
    spikes[(i-1)*1ms .<= t,i] = 1
end

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig)

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)

y,a,mi = adaptmi(identity,3spikes,
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)


# how does this look over a longer time period?

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

# what about when all inputs turn on simultaneously?

spikes = ones(length(t),20)

y,a,mi = adaptmi(identity,noise(3spikes,100ms,0.5),
                 τ_a=1.5s,c_a=3,τ_mi=50ms,c_mi=10,shape=sig,
                 W_mi=W_mi)

# what about when we use lateral, rather than global inhibition?
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

# there's some other things I need to be working on (by priority)
#
# 1. plot response of L1 to varying delta F
# 2. bistability in layer 1 (i.e. via appropriate inhibition patterns in layer 1)
# - maybe 5 goes here?
# 3. visualization of layer 2 and 3
# 4. fundamentals of bistable models, expanding to higher dimensions
# 5. figure out what the heck is up with sigmoid screwing up
#    model behavior
