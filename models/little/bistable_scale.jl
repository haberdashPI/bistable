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
dir = "../../plots/run_2018_04_16"
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))
sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

swn = drift(sweights,τ_σ = 500ms,c_σ = 0.3);
swna,m,a,x = adaptmi(swn,
                   τ_x=300ms, c_x=3.0,
                   τ_n = 2s, c_n = 5,
                   c_m=30, τ_m=350ms, W_m=scale_weighting(cs,10.0),
                   c_a=10, τ_a=3s, shape_y = x -> clamp(x,0,20))

csa = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
rplot(csa)

counts,stream1,stream2 = source_count_by_threshold(csa,cutoff=2cycoct,buildup=1s,
                                                   window=1s,delta=0.25s,
                                                   return_intermediate=true)

N = length(stream1)
df = DataFrame(x=repeat(1:N,outer=2),y=[stream1; stream2],
               stream=repeat(["stream1","stream2"],inner=N))
R"""
ggplot($df,aes(x,y,color=stream)) +
  geom_line() + scale_color_brewer(palette='Set1')
"""

lplot(counts)

crsa = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crsa,ncomponents=3,window=200ms,method=:nmf,
            delta=50ms,maxiter=100,tol=1e-3)
rplot(C)


# crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
# C = cohere(crs[0s .. 5s],ncomponents=3,window=300ms,method=:nmf,
#             delta=200ms,maxiter=100,tol=1e-3)
# rplot(C)

Ct = track(C,tc=0.5s)
rplot(Ct)

AuditoryCoherence.MultiNormalStats(reshape())
Ctp = track(C,method=:prior,tc=0.8s,prior=)

# TODO: implement simple MAP tracking with a
# prior defined by the data, and then match each trial
# to either an existing source or a new source
# on the basis of that prior.

threshold = mean(Ct) + 2std(Ct)
Ctt = copy(Ct)
Ctt[Ct .< threshold] = 0
rplot(Ctt)

# TODO: implement this for time indices
mywindow_(fn,win,step,xs) =
    [fn(xs[clamp.(i:i+win,0,end),j])
     for i in 1:step:size(xs,1), j in indices(xs,2)]

axes = tuple((axisdim(Ct,ax) for ax in [Axis{:scale},Axis{:freq}])...)
sources = squeeze(maximum(Ct .> threshold,axes),axes)
counts = sum(mywindow_(maximum,20,10,sources),2)

AuditoryModel.raster_plot(Float64.(sources))
AuditoryModel.raster_plot(Float64.(mywindow_(any,20,10,sources[:,:])))

lplot(x) = R"qplot(x=$(1:length(x)),y=$x,geom='line')"
lplot(counts)

len = (findlengths(clamp.(counts,1,2)) |> mergelengths(15)) .* (Δt(Ct)*10)


# compute 50% quantile, determine which components
# at each time slice have a value over this.


# create a plot of the auditory output....
crsC = mask(crs,component(Ct,1))
xm = audiospect(crsC)
rplot(xm)


# TODO: look at spectral masking output
# to see if we can get some statistics for the low-level
# adapt/MI

########################################
# TODO: fix below to fit with new API

# # ########################################
# # looking at noise

# xl = ab(120ms,120ms,1,150,500Hz,6) |> normpower |> amplify(-10dB)
# params = AdaptMI(c_m=1,τ_m=200ms,W_m=scale_weighting2(cs,0.75),
#                  c_a=12,τ_a=2s,shape_y = x -> max(0,x),
#                  τ_σ = 100ms, c_σ = 0.5,
#                  Δt = Δt(cs));

# cs = cortical(xl;cparams...)
# scale_noise = drift(zeros(ntimes(cs),nscales(cs)),params);

# csa = similar(cs);
# time = Axis{:time}
# csa,a,m = adaptmi(csa,params) do cs_t,t,dt_cs
#   cs[time(t)] .* (1 .+ scale_noise[t,:])
# end

# # p = rplot(csa)
# p = rplot(csa,scales = scales(csa)[1:2:end])
# R"""
# ggsave($(joinpath(dir,"4_noise_adaptmi.pdf")),$p)
# """

# crsa = cortical(csa[0s .. 20s,:,250Hz .. 1kHz];rates=(2.0.^(1:5))Hz)
# cohere = CoherenceModel(AuditoryModel.Params(cs),3,window=200ms,method=:nmf,
#                         delta=150ms,maxiter=200,tol=1e-3)
# tc = cohere(crsa)
# p = rplot(cohere,tc,crsa)
# R"""
# ggsave($(joinpath(dir,"5_cohere_noise_adaptmi.pdf")),$p)
# """


# # Thought: the noise as written probably averages out to somethign very
# # minimal across all components of a scale, I should probably have
# # per scale noise
