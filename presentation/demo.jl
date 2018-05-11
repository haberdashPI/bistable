push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "plots"
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,2,50Hz,6,:ramp) |> normpower |> amplify(-10dB)

R"""
qplot(x=$((1:length(x))./ustrip(samplerate())),y=$(Array(x)),geom='line') +
  xlab('Time (s)')
ggsave($(joinpath(dir,"timeamp.pdf")))
"""

x = ab(120ms,120ms,1,4,500Hz,6,:ramp) |> normpower |> amplify(-10dB)

y = audiospect(x)
p = rplot(y)
R"""
ggsave($(joinpath(dir,"spect.pdf")))
"""

cs = cortical(y,scales=cycoct.*round.(2.0.^linspace(-2,2,3),1))
rplot(cs,fn=abs)
R"""
ggsave($(joinpath(dir,"scales.pdf")))
"""

csr = cortical(y,scales=cycoct.*round.(2.0.^linspace(-2,2,9),1),
               rates=[Hz.*.-round.(2.0.^linspace(0,5,11),1);
                      Hz.*round.(2.0.^linspace(0,5,11),1)])

rplot(csr[:,[14,18,20],1:3:end,:],fn=abs)
R"""
ggsave($(joinpath(dir,"scale_rates.pdf")))
"""

x = ab(120ms,120ms,1,10,500Hz,6) |> normpower |> amplify(-10dB)
y = audiospect(x)
csr = cortical(y[:,400Hz .. 800Hz],scales=cycoct.*round.(2.0.^linspace(-2,2,9),1),
               rates=default_rates)


C = cohere(csr,ncomponents=3,window=350ms,method=:nmf,
           delta=100ms,maxiter=100,tol=1e-3)
rplot(C[0s .. 2s,[3,7,9],:,:])
R"""
ggsave($(joinpath(dir,"NMF.pdf")))
"""

Ct = track(C,tc=2s)
rplot(Ct[0s .. 2s,[3,7,9],:,:])
R"""
ggsave($(joinpath(dir,"tracking.pdf")))
"""

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)
y = audiospect(x)
yn = drift(y,τ_σ = 500ms,c_σ = 0.3);
yna,a,m = adaptmi(yn,
                  τ_x=300ms, c_x=3.0,
                  τ_n = 1s, c_n = 3,
                  c_m=30,τ_m=350ms,W_m=freq_weighting(y,10.0),
                  c_a=10,τ_a=3s,shape_y = x -> max(0,x))
yna_clean = similar(yna[2s .. 12s])
yna_clean .= min.(0.3,yna[2s .. 12s])
rplot(yna_clean)
R"""
ggsave($(joinpath(dir,"freq_bistable.pdf")))
"""


swn = drift(sweights,τ_σ = 500ms,c_σ = 0.3);
swna,a,m = adaptmi(swn,
                   τ_x=300ms, c_x=3.0,
                   τ_n = 1s, c_n = 15,
                   c_m=30,τ_m=350ms,W_m=scale_weighting(cs,10.0),
                   c_a=10,τ_a=3s,shape_y = x -> max(0,x))

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)

sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-2,2,9),1))
sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))


swn = drift(sweights,τ_σ = 500ms,c_σ = 0.5);
swna,a,m = adaptmi(swn,
                   τ_x=300ms, c_x=3.0,
                   τ_n = 1s, c_n = 15,
                   c_m=30,τ_m=600ms,W_m=scale_weighting(cs,10.0),
                   c_a=10,τ_a=5s,shape_y = x -> max(0,x))
csa = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)
rplot(csa[:,[3,5,7,8,9],:],fn=abs)
R"""
ggsave($(joinpath(dir,"scale_bistable.pdf")))
"""

csr = cortical(csa[:,:,400Hz .. 800Hz],rates=default_rates)
C = cohere(csr,ncomponents=3,window=500ms,method=:nmf,
           delta=250ms,maxiter=100,tol=1e-3)

Ct = track(C,tc=2s)
rplot(Ct)
R"""
ggsave($(joinpath(dir,"scale_bistable_nmf.pdf")))
"""

rplot(Ct[2.5s .. 5s,1:2:end,:,:])
R"""
ggsave($(joinpath(dir,"scale_bistable_nmf_chunk1.pdf")),width=4,height=4)
"""


rplot(Ct[6.5s .. 8.5s,1:2:end,:,:])
R"""
ggsave($(joinpath(dir,"scale_bistable_nmf_chunk2.pdf")),width=4,height=4)
"""


rplot(Ct[10.75s .. 11s,1:2:end,:,:])
R"""
ggsave($(joinpath(dir,"scale_bistable_nmf_chunk3.pdf")),width=4,height=4)
"""
