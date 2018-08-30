push!(LOAD_PATH,joinpath("packages"))
using JLD2
using AuditoryModel
using DataFrames
using AuditoryCoherence
using AxisArrays
using RCall

include("util/stim.jl")
include("util/peaks.jl")
include("util/lengths.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_$(Date(now()))"
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))

#  sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     #  axes(cs,Axis{:time}),
                     #  axes(cs,Axis{:scale}))

#  swn = drift(sweights,τ_σ = 500ms,c_σ = 0.3);
#  swna,a,m = adaptmi(swn,c_m=5,τ_m=350ms,W_m=scale_weighting(cs,1.0),
                   #  c_a=30,τ_a=3s,shape_y = x -> max(0,x),
                   #  c_x=1.0,τ_x=300ms)

#  csa = similar(cs);
#  csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
#  rplot(csa)


# @save "scratch.jld2" cs csa C
# @load "scratch.jld2" cs csa C

#  Ct = track(C,tc=0.8s)
#  rplot(Ct)

#  early = C[0s .. 4s];
#  p = reshape(mean(early,4),size(early,1),:)
#  early_prior = AuditoryCoherence.MultiNormalStats(p);
#  Ctp = track(C,method=:prior,tc=1s,prior=early_prior,thresh=1e-1)
#  rplot(Ctp)
#  alert()

crs = cortical(cs[:,:,400Hz .. 1200Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
Cna = cohere(crs,ncomponents=3,window=200ms,method=:nmf,
             delta=50ms,maxiter=100,tol=1e-3)

early = Cna[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
early_prior_iso = AuditoryCoherence.IsoMultiNormalStats(p,10);
freq_prior = AuditoryCoherence.BinomialCond(
  :old => AuditoryCoherence.Beta(1.9,2.0),
  :new => AuditoryCoherence.Beta(0.1,2.0)
);

# fused
Cna_t,sources,source_sds,lps_100 =
  track(Cna,method=:prior,tc=100ms,source_prior=early_prior_iso,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)
alert()

R"""
qplot(x=$(ustrip.(times(lps_100))),y=$(lps_100.data),geom='line')
"""

Cna_t,sources,source_sds,lps_250 =
  track(Cna,method=:prior,tc=250ms,source_prior=early_prior_iso,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)
alert()


R"""
qplot(x=$(ustrip.(times(lps_250))),y=$(lps_250.data),geom='line')
"""

df = DataFrame(x=[ustrip.(times(lps_100));ustrip.(times(lps_250))],
               y=[Array(lps_100);Array(lps_250)],
               window=[fill(100.0,ntimes(lps_100));fill(250.0,ntimes(lps_250))])

R"""
ggplot($df,aes(x,y,color=factor(window))) + geom_line() +
  scale_color_brewer(palette='Set1')
"""

# split
Cna_t,sources,source_sds =
  track(Cna,method=:prior,tc=1s,source_prior=early_prior_iso,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)
alert()

# early_prior_iso.S *= 10;

# Cna_t,sources,source_sds =
#   track(Cna,method=:prior,tc=250ms,source_prior=early_prior_iso,
#         freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
# alert()

# how does the time constant fair with a more naive prior?
p = AuditoryCoherence.IsoMultiNormalStats(10.0,10,Cna)

# fused
Cna_t,sources,source_sds,lps =
  track(Cna,method=:prior,tc=100ms,source_prior=p,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)
alert()

# split
Cna_t,sources,source_sds =
  track(Cna,method=:prior,tc=2s,source_prior=p,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)

# split
Cna_t,sources,source_sds =
  track(Cna,method=:prior,tc=500ms,source_prior=p,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)

# what about a slighlty less naive prior
μ = mean(early)
p = AuditoryCoherence.IsoMultiNormalStats(2μ,10,Cna)

Cna_t,sources,source_sds =
  track(Cna,method=:prior,tc=50ms,source_prior=p,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
rplot(Cna_t)
alert()

# OKAY, we can do this with the very naive prior,
# that seems fine, next step
# is to look at multiple tracks at once, and
# examine their logprob
# ALTHOUGH: there is something weird about
# sources are being added with these priors that
# looks different (it's grouping the off dominant tone
# and there is no component present when the dominant tone
# isn't playing)
# If we use the non naive prior, we can get fusing
# when we use a very short time constant


# ISSUE: both parameters work, SD and time constant.
# Question: are these "trivial" in the sesne that any sound
# would group? with a short enough time constant, the prior
# dominates. However, *when* this occurs likely depends
# on the nature of the stimulus, so that might be okay...

# Cna_t,sources,source_sds =
#   track(Cna,method=:prior,tc=250ms,source_prior=early_prior_iso,
#         freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
# alert()


# rplot(Cna)
# R"""
# ggsave($(joinpath(dir,"1_cohere_before_tracking.pdf")))
# """

rplot(Cna_t)
# R"""
# ggsave($(joinpath(dir,"2_cohere_after_tracking_large_SD.pdf")))
# """

# early = Cna[0s .. 4s];
# p = reshape(mean(early,4),size(early,1),:)
# early_prior_iso = AuditoryCoherence.IsoMultiNormalStats(p,10);
# freq_prior = AuditoryCoherence.BinomialCond(
#   :old => AuditoryCoherence.Beta(1.9,2.0),
#   :new => AuditoryCoherence.Beta(0.1,2.0)
# );
# early_prior_iso.S *= 10;

# Cna_t,sources,source_sds =
#   track(Cna,method=:prior,tc=250ms,source_prior=early_prior_iso,
#         freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
# alert()

# rplot(Cna_t)
# R"""
# ggsave($(joinpath(dir,"3_cohere_after_tracking_small_SD.pdf")))
# """

# TODO: the time constant, the prior variance, and the prior strength
# all influence whehther the sources divide, now let's create
# an ensemble with different values for these hyperparameters

# TODO: can I pick priors without using the actual data?
# maybe use a reasonable range? (what is the actually variance?)
alert()
