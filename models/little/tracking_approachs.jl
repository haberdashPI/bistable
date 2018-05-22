push!(LOAD_PATH,joinpath(@__DIR__,"packages"))
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
sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

swn = drift(sweights,τ_σ = 500ms,c_σ = 0.3);
swna,a,m = adaptmi(swn,c_m=5,τ_m=350ms,W_m=scale_weighting(cs,1.0),
                   c_a=30,τ_a=3s,shape_y = x -> max(0,x),
                   c_x=1.0,τ_x=300ms)

csa = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
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
early_prior_iso = AuditoryCoherence.IsoMultiNormalStats(p);
Cna_t = track(Cna,method=:prior,tc=1s,prior=early_prior_iso,thresh=1e-1)
rplot(Cna_t)
alert()

crsa = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crsa,ncomponents=3,window=200ms,method=:nmf,
            delta=50ms,maxiter=100,tol=1e-3)
#  rplot(C)

early = C[0s .. 4s];
p = reshape(mean(early,4),size(early,1),:)
early_prior_iso = AuditoryCoherence.IsoMultiNormalStats(p);
Ct = track(C,method=:prior,tc=1s,prior=early_prior_iso,thresh=1e-1)
rplot(Ct)
alert()
# last block above leads to this error
#=

ERROR: LoadError: MethodError: no method matching AuditoryCoherence.IsoMultiNormalStats{Float64}(::Arr
Stacktrace:
 [1] AuditoryCoherence.IsoMultiNormalStats(::Array{Float64,2}, ::Int64, ::Int64, ::Int64) at /Users/da
jl:91
 [2] AuditoryCoherence.IsoMultiNormalStats(::Array{Float64,2}) at /Users/davidlittle/MEGA/Research/sof
 [3] include_from_node1(::String) at ./loading.jl:576
 [4] include(::String) at ./sysimg.jl:14
while loading /Users/davidlittle/MEGA/Research/software/bistable/models/little/tracking_approachs.jl,

=#
