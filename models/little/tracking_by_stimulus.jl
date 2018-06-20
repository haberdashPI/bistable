push!(LOAD_PATH,joinpath("packages"))
using JLD2
using AuditoryModel
using DataFrames
using AuditoryCoherence
using AxisArrays
using RCall

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_$(Date(now()))"
isdir(dir) || mkdir(dir)

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

Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms,500ms],
               thresh=1e-1,source_priors = priors,
               max_sources = 3, freq_prior = freq_prior)
alert()

p1 = rplot(Ct[1])
p2 = rplot(Ct[2])

R"""
library(cowplot)
p = plot_grid($p1 + ggtitle(paste("Time constant =",$(string(tcs[1])))),
              $p2 + ggtitle(paste("Time constant =",$(string(tcs[2])))),
              nrow=2)
save_plot($(joinpath(dir,"1_stream_by_time_constant.pdf")),p,
          base_aspect_ratio=1.3,nrow=2,ncol=1)
"""

df1 = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(times(lp)),lps)...),
               tc = repeat(ustrip.(tcs),inner=length(lps[1])),
               prior = repeat(ps,inner=length(lps[1])),
               delta_f = 6)

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
for c in [1] #[0.5,1,2]
  cur_prior = deepcopy(base_prior)
  cur_prior.S *= c;
  push!(priors,cur_prior)
end

Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms,500ms],
                      thresh=1e-1,source_priors = priors,
                      max_sources = 3, freq_prior
                      = freq_prior)

df2 = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(times(lp)),lps)...),
               tc = repeat(ustrip.(tcs),inner=length(lps[1])),
               prior = repeat(ps,inner=length(lps[1])),
               delta_f = 12)


x = ab(120ms,120ms,1,25,500Hz,2) |> normpower |> amplify(-10dB)

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

Ct,lps,tcs,ps = track(C,method=:multi_prior,tcs = [100ms,250ms,500ms],
                      thresh=1e-1,source_priors = priors,
                      max_sources = 3, freq_prior = freq_prior)

df3 = DataFrame(logpdf = vcat(Array.(lps)...),
               time = vcat(map(lp -> ustrip.(times(lp)),lps)...),
               tc = repeat(ustrip.(tcs),inner=length(lps[1])),
               prior = repeat(ps,inner=length(lps[1])),
               delta_f = 2)


df = [df1; df2; df3]

R"""
dfp = subset($df,tc %in% c(0.1,0.25))
p = ggplot(dfp,aes(x=time,y=logpdf,color=factor(tc))) +
  facet_wrap(~paste("Delta f =",delta_f)) +
  scale_color_brewer(palette='Set1') + geom_line()
save_plot($(joinpath(dir,"2_object_logpdf_by_delta_f.pdf")),p,
  base_aspect_ratio=1.1,nrow=1,ncol=3)
"""

