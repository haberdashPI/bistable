push!(LOAD_PATH,"packages")
using AuditoryModel
using DataFrames
using AuditoryCoherence
using AxisArrays
using RCall
using JLD

include("util/stim.jl")
include("util/peaks.jl")
include("util/lengths.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_04_16"
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))
sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

swn = drift(sweights,τ_σ = 500ms,c_σ = 0.3);
swna,a,m = adaptmi(swn,
                   c_m=5,τ_m=350ms,W_m=scale_weighting2(cs,1.0),
                   c_e=1.5,τ_e = 300ms,
                   c_a=30,τ_a=3s,shape_y = x -> max(0,x),α=3.0)

csa = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
rplot(csa)

#########################################
# find object count by peak picking

# first, a little experiment on a single slice
lplot(x) = R"""
qplot(x=1:$(length(x)),y=$(vec(x)),geom='line')
"""

slice1 = vec(abs.(mean(Array(csa[0s .. 1s,0.5cycoct,:]),1)))

R"""
qplot(x=1:$(length(slice1)),y=$slice1,geom='line') +
  geom_vline(xintercept=$(find_peaks(slice1)))
"""

slice2 = vec(abs.(mean(Array(csa[5s .. 6s,3.1cycoct,:]),1)))

R"""
qplot(x=1:$(length(slice2)),y=$slice2,geom='line') +
  geom_vline(xintercept=$(find_peaks(slice2)))
"""

win = abs.(squeeze(mean(Array(csa[0s .. 1s,:,:]),1),1))
guess_source_count(win)

win = abs.(squeeze(mean(Array(csa[5s .. 6s,:,:]),1),1))
guess_source_count(win)

counts = source_count_by_peaks(csa)

# TODO: things seem to be make sense for high cyc/oct, what about low?
p1 = rplot(source_bumps(csa[:,2.4cycoct .. 3.1cycoct,:]));
p2 = rplot(source_peaks(csa[:,2.4cycoct .. 3.1cycoct,:],stage=:final));

R"""
plot_grid($p1,$p2,nrow=2)
"""

p1 = rplot(source_bumps(csa[:,0.5cycoct .. 0.6cycoct,:]));
p2 = rplot(source_peaks(csa[:,0.5cycoct .. 0.6cycoct,:],stage=:final));

R"""
plot_grid($p1,$p2,nrow=2)
"""
counts = source_count_by_peaks(csa)
R"""
qplot(x=$(ustrip.(times(counts))),y=$(Array(counts)),geom='line')
"""

lens,vals = findlengths(counts)

########################################
# find object count by scales

# find a good threshold
threshold = mean(abs.(csa))
csat = copy(csa)
csat[abs.(csa) .< threshold] = 0
rplot(csat)

# okay... seems like a good threshold
two_stream_scale = 2cycoct

onestream = sum(max.(0,abs.(csa[:,0cycoct .. two_stream_scale,:]) .- threshold),
                [axisdim(csa,Axis{:scale}) axisdim(csa,Axis{:freq})])
twostream = sum(max.(0,abs.(csa[:,two_stream_scale .. 4cycoct,:]) .- threshold),
                [axisdim(csa,Axis{:scale}) axisdim(csa,Axis{:freq})])

window(fn,win,xs::Array) = [fn(xs[clamp.(i:i+win,0,end)])
                            for i in 1:win:length(xs)]

onemean = window(mean,20,onestream)
twomean = window(mean,20,twostream)
N = length(onemean)
ts = times(csa)[1:20:ntimes(csa)]

onemean ./= maximum(onemean[ts .> 1s])
twomean ./= maximum(twomean)

df = DataFrame(time =  ustrip.([ts; ts]),
               sum = [vec(onemean);vec(twomean)],
               kind = [fill("onestream",N);fill("twostream",N)])
R"""
ggplot($df,aes(x=time,y=sum,color=kind)) + geom_line() +
    scale_color_brewer(palette='Set1')
"""

findlengths(x) = diff([0; find(diff(x) .!= 0) .+ 1; length(x)])
lengths = findlengths(onemean .> twomean) .* (ts[2] - ts[1])

crsa = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crsa,ncomponents=3,window=200ms,method=:nmf,
            delta=50ms,maxiter=100,tol=1e-3)
rplot(C)


# crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
# C = cohere(crs[0s .. 5s],ncomponents=3,window=300ms,method=:nmf,
#             delta=200ms,maxiter=100,tol=1e-3)
# rplot(C)

Ct = track(C,tc=0.8s)
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
