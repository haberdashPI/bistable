push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using DataFrames
using AxisArrays
using RCall

include("util/stim.jl")
include("util/lengths.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
lplot(x) = R"qplot(x=1:$(length(x)),y=$(Array(x)),geom='line')"

x        = ab(120ms,120ms,1,40,500Hz,6) |> normpower |> amplify(-10dB)
sp       = audiospect(x)
cs       = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-2,2.5,9),1))

sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

bound(x,min,max) = 1/(1+exp(-4((x - min)/(max - min) - 0.5)))
saturated  = AxisArray(bound.(sweights,0.0,0.05),axes(sweights)...)
swn        = drift(saturated,τ_σ = 500ms,c_σ = 0.3);
swna,m,a,x = adaptmi(
  swn, τ_x=500ms, c_x=3.0, τ_n=2s, c_n=5,
  c_m=30, τ_m=350ms, W_m=scale_weighting(cs,15,6),
  c_a=10, τ_a=50s, shape_y=x->clamp(x,0,20)
)
# rplot(swna)
# alert()

csa  = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)

crs = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
Ca = cohere(crs,ncomponents=3,window=150ms,method=:nmf,skipframes=1,
            delta=75ms,maxiter=100,tol=1e-3)
# rplot(Ca)
# alert()

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=150ms,method=:nmf,skipframes=2,
           delta=75ms,maxiter=100,tol=1e-3)
# rplot(C)
# alert()

# TODO: overall - get tracking to not always sum, when sum is possible
# for the bistable scales

# still seem to be problems in splitting high scales and
# fusing small ones at the same time for the bistable scales
#
# one possibility is that we could have an assumed ridge in the
# correlation matrix, that would make the broader scale
# components appear more similar, this might allow for gradual
# shifts in the distribution across channels
#
# that's what I'm working on now

# THOUGHT: might be able to handle larger scales in the ridge by excluding some
# number of indices for the lower scales (so we can "expand" the ridge at these
# lower scales without increasing computational complexity)
#
# NOTE: right now the ridge doesn't seem necessary, though
# the formulation defined by the ridge code is working better
# probably because of a math error in isonorm
#
# TODO: right now I'm trying to get a better windowing
# approach, that can depend less on older samples

dist(a,b) = (a[1] - b[1])^2 / (0.5^2) + (a[2] - b[2])^2 / (0.5^2)

Ct,source,sourceS,lp,tracks = track(
  C,
  method=:prior,
  tc=1.5s,
  source_prior=ridgenorm(C[0s .. 4s],10,scale=0.25,freq=0.25),
  freq_prior=freqprior(0,2),
  max_sources=4,
  unmodeled_prior=0
)
rplot(Ct)

# great, with a slightly larger window, that now seems to work

# next steps:
# 1. try ridge prior on all the below
# code, confirm that we're in the same place as before
# (since with out the ridge it should be equivalent)

# fusing at low scales
Cw = C[0s .. 4s,1:3,:,:]
ridgep = ridgenorm(C[0s .. 2s,1:3,:,:],10,
                 scale=0.25,freq=0.25)
ridgep.S .*= 1
# ridgep.μ .= 0
Ct,source,sourceS,lp,tracks = track(Cw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# splitting at high scales
Cw = C[0s .. 4s,6:9,:,:]
ridgep = ridgenorm(C[0s .. 2s,6:9,:,:],10,
                 scale=0.25,freq=0.25)
Ct,source,sourceS,lp,tracks = track(Cw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of fusing for low scale (with adaptmi)
Caw = Ca[0s .. 3s,1:3,:,:]
ridgep = ridgenorm(C[0s .. 1s,1:3,:,:],10,scale=0.25,freq=0.25, thresh=0.3)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of splitting at high scales (with adaptmi)
# not working very well....
Caw = Ca[6s .. 10s,7:9,:,:]
ridgep = ridgenorm(C[0s .. 1s,7:9,:,:],10, scale=0.25,freq=0.25, thresh=0.3)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of splitting at high scales (with adaptmi)
# not working very well....
Caw = Ca[9s .. 13s,8:9,:,:]
ridgep = ridgenorm(C[0s .. 1s,8:9,:,:],10, scale=0.25,freq=0.25,thresh=0.1)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# ugh... the above isn't working super well,
# it looks like it's fragile to differences in how
# things start because with the stimulus I was using
# 9s .. 11s doesn't work, but 6s .. 10s does

# TODO: this is where I left off
# think about the next steps
# (I think I need to do more than
# just change windowing of individual sources
# there is still a lot of reliance on prior state
# that I don't like, probalby about
# the greedy approach I'm using).

# okay, we're in a pretty good place now
# let's try looking across all scales, and then all times

# example of splitting (with adaptmi)
# not working very well at first, but....

# NOTE: the thresh needs to be high enough
# to make sure the all the zeros in non-responses scales
# don't outweigh what's happening with the stronger scales

Caw = Ca[7.5s .. 12s,1:9,:,:]
ridgep = ridgenorm(C[0s .. 1s,1:9,:,:],10, scale=0.25,freq=0.25,thresh=0.3)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of fusing (with adaptmi)
Caw = Ca[0s .. 3s,:,:,:]
ridgep = ridgenorm(C[0s .. 1s,:,:,:],10,scale=0.25,freq=0.25,thresh=0.3)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# can we handle this fragility with start window
# by have occasional restarts, where we check the probability
# of that restart? then just reset when the restart is better

# to test this, look at lp for both windows

# example of entire Ca
Caw = Ca[6.5s .. 12s]
ridgep = ridgenorm(C[0s .. 1s,:,:,:],10,scale=0.25,freq=0.25,thresh=0.3)
Ct,source,sourceS,lp1,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of entire Ca
Caw = Ca[7.5s .. 12s]
ridgep = ridgenorm(C[0s .. 1s,:,:,:],10,scale=0.25,freq=0.25,thresh=0.3)
Ct,source,sourceS,lp2,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

df = DataFrame(x = [ustrip.(times(lp1)); ustrip.(times(lp2))],
               y = [lp1; lp2],
               window = [fill("from6.5",length(lp1));
                         fill("from7.5",length(lp2))])
R"""
ggplot($df,aes(x,y,color=window)) + geom_line() +
  scale_color_brewer(palette='Set1')
"""

# STOPPED HERE
################################################################################

# doesn't work so well now... is this because of the transitions?
# or a longer effect on the prior?

# example of entire Ca
Caw = Ca[5.5s .. 7.5s]
ridgep = ridgenorm(C[0s .. 1s,:,:,:],10,scale=0.25,freq=0.25,thresh=0.3)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# seems to be due to the transition. does this happen with
# just the specific, active scales?

# example of entire Ca
Caw = Ca[5.5s .. 10s,6:9,:,:]
ridgep = ridgenorm(C[0s .. 1s,6:9,:,:],10,scale=0.25,freq=0.25,thresh=0.3)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# okay, so yeah, it works a little better, but the bottom line seems to be that
# there is too much dependences on older state, because just 500ms added to the
# front really screws things up, what we need is a sharper
# dropoff than the obvious decay process provides. I could do windowing or I
# could change the decay formula

################################################################################
