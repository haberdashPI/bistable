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
cs       = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))
sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

# TODO: think through these parameters
# NOTE: one very simple thing to do would be to saturate
# weights, so that the lower and higher scales have similar
# input to the adapt/MI process

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

crs = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=150ms,method=:nmf,skipframes=2,
            delta=100ms,maxiter=100,tol=1e-3)
rplot(C)
alert()

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=100ms,method=:nmf,skipframes=2,
           delta=50ms,maxiter=100,tol=1e-3)
rplot(C)
alert()

Ct,source,sourceS,lp,tracks = track(C,method=:prior,tc=2s,
                                    source_prior=isoprior(0.05,10),
                                    freq_prior=freqprior(0.9,2),thresh=1e-1,
                                    max_sources=5)

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
#             source_prior=isoprior(0.1,10),
#             freq_prior=freqprior(0,2),thresh=1e-1,max_sources=5)
# rplot(Ct1)

Ct,source,sourceS,lp,tracks = track(C[:,:,:,1:3],method=:prior,tc=1s,
                              source_prior=isoprior(0.1,2),
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
                              source_prior=isoprior(0.01,10),
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
                              source_prior=isoprior(0.01,10),
                              freq_prior=freqprior(0,2),thresh=1e-1,
                              max_sources=5)

Cnt = track(C2n,method=:simple,tc=1s)
