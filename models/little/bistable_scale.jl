push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using ProgressMeter
using DataFrames
using AxisArrays
using RCall
using DSP

include("util/stim.jl")
include("util/lengths.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
lplot(x) = R"qplot(x=1:$(length(x)),y=$(Array(x)),geom='line')"

x        = ab(120ms,120ms,1,40,500Hz,6) |> normpower |> amplify(-10dB)
sp       = audiospect(x)
cs       = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,3,8),1),
                   bandonly=true)

sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

bound(x,min,max) = 1/(1+exp(-4((x - min)/(max - min) - 0.5)))
saturated  = AxisArray(bound.(sweights,0.0,0.05),axes(sweights)...)
swn        = drift(saturated,τ_σ = 500ms,c_σ = 0.3);
swna,m,a,x = adaptmi(
  swn, τ_x=500ms, c_x=3.0,
  c_m=30, τ_m=350ms, W_m=scale_weighting(cs,15,6),
  c_a=10, τ_a=50s, shape_y=x->clamp(x,0,20)
)
low = digitalfilter(Lowpass(1.5;fs=ustrip(uconvert(s,1/Δt(swna)))),Butterworth(3))
swna_low = AxisArray(max.(0,filt(low,swna)),Axis{:time}(times(swna)))
rplot(swna_low)
# rplot(swna)
# alert()

csa  = similar(cs);
csa .= sqrt.(abs.(cs) .* swna_low) .* exp.(angle.(cs)*im)

# without adaptmi
crs = cortical(cs[:,:,400Hz .. 1200Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz],
               bandonly=true)
C = cohere(crs[0s .. 2s],ncomponents=3,window=100ms,method=:nmf,skipframes=0,
           delta=75ms,maxiter=100,tol=1e-3)


# with adaptmi
crs = cortical(csa[:,:,400Hz .. 1200Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz],
              bandonly=true)
Ca = cohere(crs,ncomponents=3,window=100ms,method=:nmf,skipframes=0,
           delta=75ms,maxiter=100,tol=1e-3)

ridgep = ridgenorm(C[0s .. 1s,:,:,:],10,scale=0.25,freq=0.25,thresh=0.05)
Ct,lp,source,sourceS,tracks = track(Ca,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

strs = map_windowing(Ct,length=500ms,step=250ms) do window
  component_means(window)
end
strmat = AxisArray(hcat(strs...)',Axis{:time}(times(strs)))
quartz(); rplot(strmat)

progress = Progress(length(windowing(Ct,length=500ms,step=250ms)))
ratios = map_windowing(Ct,length=500ms,step=250ms) do window
  strengths = sort(component_means(window),rev=true)
  ProgressMeter.next!(progress)
  strengths[1] / sum(strengths[2:end])
end
p = rplot(ratios)
R"""
$p + geom_hline(yintercept=2,linetype=2)
"""

percept_lengths(AxisArray(ratios .> 2,axes(ratios,1)),1s)

# okay, that works, let's try this with the multi-tracking method

# first verify that different time constants still yield different results
# NOTE: that's not working, but I think I can go back to using the prior
# because the method will be applied to all conditions

Cprior = C[0s .. 1s,:,:,:]
Cprior .*= 10
ridgep = ridgenorm(Cprior,10,scale=0.25,freq=0.25,thresh=0.05)

ridgep = ridgenorm(0.05,4,size(C,2,3),rdist(scale=0.25,freq=0.25))
Ct,lp,source,sourceS,tracks = track(Ca,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

ridgep = ridgenorm(0.5,,4,size(C,2,3),rdist(scale=0.25,freq=0.25))
Ct,lp,source,sourceS,tracks = track(Ca,method=:prior,tc=2s,
                                    source_prior=ridgep,
                                    freq_prior=freqprior(0,2),
                                    max_sources=4,unmodeled_prior=0)
rplot(Ct)

# now try with the multi prior

sds = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]
ridgeps = [ridgenorm(x,4,size(C,2,3),scale=0.25,freq=0.25) for x in sds]
tracks = track(Ca,method=:multi_prior,tcs=[2s],
               source_priors=ridgeps,freq_prior=freqprior(0,2),
               max_sources=4,unmodeled_prior=0)

lp_index = 4 # (will be 2)
df = DataFrame(sd = repeat(sds,inner=ntimes(tracks[1][1])),
               lp = vcat(map(t -> t[4],tracks)...),
               time = vcat(map(t -> ustrip.(uconvert.(s,times(t[1]))),tracks)...))
R"""
library(ggplot2)
ggplot($df,aes(x=time,y=lp,color=factor(sd))) +
  geom_line() + scale_color_brewer(palette='Reds')
"""


