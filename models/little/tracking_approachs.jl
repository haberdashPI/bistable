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

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
Cna = cohere(crs,ncomponents=3,window=200ms,method=:nmf,
             delta=50ms,maxiter=100,tol=1e-3)

early = Cna[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
early_prior_iso = AuditoryCoherence.IsoMultiNormalStats(p,10);
freq_prior = AuditoryCoherence.BinomialCond(
  :old => AuditoryCoherence.Beta(1.9,2.0),
  :new => AuditoryCoherence.Beta(0.1,2.0)
);
early_prior_iso.S *= 1;

Cna_t,sources,source_sds =
  track(Cna,method=:prior,tc=250ms,source_prior=early_prior_iso,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
alert()

rplot(Cna)
R"""
ggsave($(joinpath(dir,"1_cohere_before_tracking.pdf")))
"""

rplot(Cna_t)
R"""
ggsave($(joinpath(dir,"2_cohere_after_tracking_large_SD.pdf")))
"""

early = Cna[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
early_prior_iso = AuditoryCoherence.IsoMultiNormalStats(p,10);
freq_prior = AuditoryCoherence.BinomialCond(
  :old => AuditoryCoherence.Beta(1.9,2.0),
  :new => AuditoryCoherence.Beta(0.1,2.0)
);
early_prior_iso.S *= 10;

Cna_t,sources,source_sds =
  track(Cna,method=:prior,tc=250ms,source_prior=early_prior_iso,
        freq_prior=freq_prior,thresh=1e-1,max_sources = 3)
alert()

rplot(Cna_t)
R"""
ggsave($(joinpath(dir,"3_cohere_after_tracking_small_SD.pdf")))
"""

# TODO: the time constant, the prior variance, and the prior strength
# all influence whehther the sources divide, now let's create
# an ensemble with different values for these hyperparameters

# TODO: can I pick priors without using the actual data?
# maybe use a reasonable range? (what is the actually variance?)
alert()
