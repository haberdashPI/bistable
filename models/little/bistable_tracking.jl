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
for c in [1] #[0.5,1,2]
  cur_prior = deepcopy(base_prior)
  cur_prior.S *= c;
  push!(priors,cur_prior)
end

Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms],
               thresh=1e-1,source_priors = priors,
               max_sources = 3, freq_prior = freq_prior)

df = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(uconvert.(s,times(lp))),lps)...),
               tc = repeat(ustrip.(uconert.(s,tcs)),inner=length(lps[1])),
               prior = repeat(ps,inner=length(lps[1])))

R"""
ggplot($df,aes(x=time,y=logpdf,color=factor(1000*tc))) + geom_line() +
  facet_wrap(~paste("prior =",round(prior,3)))
"""

oweights = AxisArray(hcat(lps...),axes(C,Axis{:time}),Axis{:prior}(1:length(lps)))
oweights .-= minimum(oweights[0s .. 4s])
oweights ./= 1000

function even_weightings(n)
  W = 2(1/(n-1))*(ones(n,n) - I)
  x -> W*x
end

own = drift(oweights,τ_σ = 500ms,c_σ = 0.1);
owna,m,a,x,n = adaptmi(own,
                   τ_x=300ms, c_x=1.5,
                   τ_n = 1s, c_n = 5,
                   c_m=20, τ_m=350ms, W_m=even_weightings(length(lps)),
                   c_a=10, τ_a=3s, shape_y = x -> clamp(x,0,20))

rplot(owna)
rplot(m)
rplot(a)
rplot(x)
rplot(n)
