push!(LOAD_PATH,joinpath("packages"))
using JLD2
using AuditoryModel
using DataFrames
using AuditoryCoherence
using AxisArrays
using RCall

quartz() = R"quartz()"

include("util/stim.jl")
include("util/peaks.jl")
include("util/lengths.jl")

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=200ms,method=:nmf,
             delta=50ms,maxiter=100,tol=1e-3)

early = C[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
base_prior = AuditoryCoherence.IsoMultiNormalStats(p,10);
freq_prior = AuditoryCoherence.BinomialCond(
  :old => AuditoryCoherence.Beta(1.9,2.0),
  :new => AuditoryCoherence.Beta(0.1,2.0)
);

priors = []
for c in [0.5,1,2]
  cur_prior = deepcopy(base_prior)
  cur_prior.S *= c;
  push!(priors,cur_prior)
end

Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms,500ms],
               thresh=1e-1,source_priors = priors,
               max_sources = 3, freq_prior = freq_prior)

df = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(times(lp)),lps)...),
               tc = repeat(ustrip.(tcs),inner=length(lps[1])),
               prior = repeat(ps,inner=length(lps[1])))

R"""
ggplot($df,aes(x=time,y=logpdf,color=factor(1000*tc))) + geom_line() +
  facet_wrap(~paste("prior =",round(prior,3)))
"""

oweights = AxisArray(hcat(lps...),axes(C,Axis{:time}),Axis{:prior}(1:length(lps)))
oweights .-= minimum(oweights[0s .. 4s])


################################################################################
# how does these same graphs look for stimuli that should consistently stream?

x = ab(120ms,120ms,1,25,500Hz,12) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=200ms,method=:nmf,
             delta=50ms,maxiter=100,tol=1e-3)

early = C[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
base_prior = AuditoryCoherence.IsoMultiNormalStats(p,10);
freq_prior = AuditoryCoherence.BinomialCond(
  :old => AuditoryCoherence.Beta(1.9,2.0),
  :new => AuditoryCoherence.Beta(0.1,2.0)
);

priors = []
for c in [0.5,1,2]
  cur_prior = deepcopy(base_prior)
  cur_prior.S *= c;
  push!(priors,cur_prior)
end


Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms,500ms],
                      thresh=1e-1,source_priors = priors,
                      max_sources = 3, freq_prior
                      = freq_prior)

df = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(times(lp)),lps)...),
               tc =
               repeat(ustrip.(tcs),inner=length(lps[1])),
               prior =
               repeat(ps,inner=length(lps[1])))

quartz()

R"""
ggplot($df,aes(x=time,y=logpdf,color=factor(1000*tc))) + geom_line() +
  facet_wrap(~paste("prior =",round(prior,3)))
"""

# what about when the sounds should 'always' fuse?
x = ab(120ms,120ms,1,25,500Hz,12) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=200ms,method=:nmf,
             delta=50ms,maxiter=100,tol=1e-3)

early = C[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
base_prior = AuditoryCoherence.IsoMultiNormalStats(p,10);
freq_prior = AuditoryCoherence.BinomialCond(
  :old => AuditoryCoherence.Beta(1.9,2.0),
  :new => AuditoryCoherence.Beta(0.1,2.0)
);

priors = []
for c in [0.5,1,2]
  cur_prior = deepcopy(base_prior)
  cur_prior.S *= c;
  push!(priors,cur_prior)
end


Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms,500ms],
                      thresh=1e-1,source_priors = priors,
                      max_sources = 3, freq_prior
                      = freq_prior)

df = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(times(lp)),lps)...),
               tc =
               repeat(ustrip.(tcs),inner=length(lps[1])),
               prior =
               repeat(ps,inner=length(lps[1])))

quartz()

R"""
ggplot($df,aes(x=time,y=logpdf,color=factor(1000*tc))) + geom_line() +
  facet_wrap(~paste("prior =",round(prior,3)))
"""
