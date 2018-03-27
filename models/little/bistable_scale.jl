push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using AxisArrays
using RCall
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_03_10"
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-20dB)

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

csa = similar(cs)
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
rplot(csa)

# TODO: verify that the new simple tracking implementation
# works for the more basic, non adaptmi input.

crsa = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crsa[0s .. 10s],ncomponents=3,window=300ms,method=:nmf,
            delta=200ms,maxiter=100,tol=1e-3)
rplot(C)

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs[0s .. 5s],ncomponents=3,window=300ms,method=:nmf,
            delta=200ms,maxiter=100,tol=1e-3)
rplot(C)

Ct = track(C,tc=0.6s)
rplot(Ct)

crsC = mask(crs,component(C,1))
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
