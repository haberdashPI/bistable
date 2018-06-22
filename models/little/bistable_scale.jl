push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using DataFrames
using AxisArrays
using RCall

include("util/stim.jl")
include("util/peaks.jl")
include("util/lengths.jl")
include("util/threshold.jl")

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
saturated  = AxisArray(bound.(sweights,0.0,0.08),axes(sweights)...)
swn        = drift(saturated,τ_σ = 500ms,c_σ = 0.3);
swna,m,a,x = adaptmi(
  swn, τ_x=100ms, c_x=3.0, τ_n=2s, c_n=5,
  c_m=30, τ_m=350ms, W_m=scale_weighting(cs,15,6),
  c_a=10, τ_a=50s, shape_y=x->clamp(x,0,20)
)

# rplot(swna)
# alert()

csa  = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)

# crs = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
# Ca = cohere(crs,ncomponents=3,window=150ms,method=:nmf,skipframes=1,
#             delta=100ms,maxiter=100,tol=1e-3)

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=150ms,method=:nmf,skipframes=2,
           delta=100ms,maxiter=100,tol=1e-3)
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

# TODO: testing new ridge normal prior currently doing some sanity checks using
# test_ridge.jl, maybe do a few more (even though it says "okay..." in the file
# ) and then try out the code below


# TODO: I should not need the threshold parameter if I'm doing
# the prior right. I can get rid of that implementation detail.

dist(a,b) = (a[1] - b[1])^2 / (0.5^2) + (a[2] - b[2])^2 / (0.5^2)

Ct,source,sourceS,lp,tracks = track(
  C,
  method=:prior,
  tc=1s,
  source_prior=ridgenorm(C[0s .. 4s],10,scale=0.25,freq=0.25),
  freq_prior=freqprior(0,2),
  thresh=1e-3,
  max_sources=4,
  unmodeled_prior=0
)

alert()
rplot(Ct)

################################################################################
# OLD STUFF (this is what convinced me to try the ridge prior)

Ct,source,sourceS,lp,tracks = track(C,method=:prior,tc=1s,
                                    source_prior=isonorm(C[0s .. 4s],10),
                                    freq_prior=freqprior(0,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
alert()
rplot(Ct)

# fusing at low scales
Cw = C[0s .. 4s,1:3,:,:]
isop = isonorm(C[0s .. 2s,1:3,:,:],1)
isop.S .*= 1
# isop.μ .= 0
Ct,source,sourceS,lp,tracks = track(Cw,method=:prior,tc=0.5s,
                                    source_prior=isop,
                                    freq_prior=freqprior(0,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# splitting at high scales
Cw = C[0s .. 4s,6:9,:,:]
isop = isonorm(C[0s .. 2s,6:9,:,:],1)
Ct,source,sourceS,lp,tracks = track(Cw,method=:prior,tc=0.5s,
                                    source_prior=isop,
                                    freq_prior=freqprior(0,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of fusing for low scale
Caw = Ca[0s .. 2s,1:3,:,:]
isop = isonorm(C[0s .. 1s,1:3,:,:],1)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=0.5s,
                                    source_prior=isop,
                                    freq_prior=freqprior(0,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of splitting at high scales
Caw = Ca[13s .. 20s,9:9,:,:]
isop = isonorm(C[0s .. 1s,9:9,:,:],1)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=0.5s,
                                    source_prior=isop,
                                    freq_prior=freqprior(0,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)


################################################################################

# example of fusing across many scales
Caw = Ca[1s .. 5s,1:3,:,:]
Caw ./= clamp.(mean(Caw,3),1e-1,Inf)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=2s,
                                    source_prior=isonorm(0.4,0.25),
                                    freq_prior=freqprior(0.9,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# example of splitting across many scales
Caw = Ca[13s .. 17s,7:9,:,:]
# Caw ./= clamp.(mean(Caw,3),1e-1,Inf)
Ct,source,sourceS,lp,tracks = track(Caw,method=:prior,tc=1s,
                                    source_prior=isonorm(0.15,0.25),
                                    freq_prior=freqprior(0.9,2),thresh=1e-3,
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

########################################
Ct,source,sourceS,lp,tracks = track(Ca[14s .. 16s],method=:prior,tc=2s,
                                    source_prior=isonorm(0.15,1),
                                    freq_prior=freqprior(0,2),thresh=1e-3,
                                    max_sources=5,unmodeled_prior=0)
alert()
rplot(Ct)

# Ct,source,sourceS,lp,tracks = track(C,method=:prior,tc=2s,
#                                     source_prior=isonorm(0.05,10),
#                                     freq_prior=freqprior(0.9,2),thresh=1e-1,
#                                     max_sources=5)

# NOTE: this analysis works quite quickly, so the
# issue with it not working with the bistable stimulus
# must be something about the more dynamic responses
# I could model some sort of delta to handle this

rplot(Ct)
alert()

# NEXT PLAN: create a continuum from
# a stimulus that works to one that doesn't
# and identify what about the stimulus is a problem

# Ct1, = track(C1,method=:prior,tc=250ms,
#             source_prior=isonorm(0.1,10),
#             freq_prior=freqprior(0,2),thresh=1e-1,max_sources=5)
# rplot(Ct1)

Ct,source,sourceS,lp,tracks = track(C[:,:,:,1:3],method=:prior,tc=1s,
                              source_prior=isonorm(0.1,2),
                              freq_prior=freqprior(0,100),thresh=0,
                              max_sources=5)
rplot(Ct)
alert()

#=
Okay... this isn't working. I need to think about either how to make a simple
change that might help or find a better way to diagnose what is going wrong.

more thoughts - It seems like we're forming sources from the
"background" source, that is the source that isn't playing
right now. What if I normalize the values before
tracking, that might help the tracking equate foreground and
background, and maybe make it easier to split
=#
####

C2n = deepcopy(C2)
C2n ./= max.(1e-2,mean(C2,[2,3]))

Cnt,source,sourceS,lp = track(C2n[12s .. 14s,9:9,:],method=:prior,tc=1s,
                              source_prior=isonorm(0.01,10),
                              freq_prior=freqprior(0,2),thresh=1e-1,
                              max_sources=5)
rplot(Cnt)

#=
STOPPED HERE!!

okay... that seemed to lead to fewer sources. Now we just have the
one merged source and noise

is there some way to visualize why other interpretations
don't win?

my first thought was to consider showing the second best tracked sources, but
this is potentially quite non-sensical, since it may vary substantially frame to
frame and the naive visualization (using the same method to build a source model
as the top interpretation) would be a function of sums across frames
=#

Cnt,source,sourceS,lp = track(C2n,method=:prior,tc=250ms,
                              source_prior=isonorm(0.01,10),
                              freq_prior=freqprior(0,2),thresh=1e-1,
                              max_sources=5)

Cnt = track(C2n,method=:simple,tc=1s)
