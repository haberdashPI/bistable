push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using AxisArrays
using RCall
using DataFrames

include("util/stim.jl")
R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_$(Date(now()))"
isdir(dir) || mkdir(dir)
lplot(x) = R"qplot(x=1:$(length(x)),y=$(vec(x)),geom='line')"

x = ab(120ms,0ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)

sp = audiospect(x)

spn = drift(sp,τ_σ = 500ms,c_σ = 0.3);

sigma = 1/2
lplot(freq_weighting(spn,sigma)(sp[1,:]))

spna,a,m,e,d = adaptmi(spn,
                       c_m=50,τ_m=150ms,W_m=freq_weighting(spn,sigma),
                       c_d=11,τ_d = 300ms,
                       c_e=10,τ_e = 200ms,
                       c_a=3,τ_a=3s,shape_y = x -> max(0,x),α=3.0)

# rplot(m[2s .. 12s])
rplot(spna[2s .. 12s])

alert()
