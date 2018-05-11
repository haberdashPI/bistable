
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
dir = "../../plots/run_" * string(Date(now()))
isdir(dir) || mkdir(dir)

x = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)
sp = audiospect(x)
cs = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-1,2,9),1))
sweights = AxisArray(squeeze(mean(abs.(cs),axisdim(cs,Axis{:freq})),3),
                     axes(cs,Axis{:time}),
                     axes(cs,Axis{:scale}))

swn = drift(sweights,τ_σ = 500ms,c_σ = 0.3);
swna,a,m = adaptmi(swn,
                   τ_x=300ms, c_x=3.0,
                   τ_n = 1s, c_n = 15,
                   c_m=30,τ_m=350ms,W_m=scale_weighting(cs,10.0),
                   c_a=10,τ_a=3s,shape_y = x -> max(0,x))

csa = similar(cs);
csa .= sqrt.(abs.(cs) .* swna) .* exp.(angle.(cs)*im)

crsa = cortical(csa[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crsa,ncomponents=3,window=200ms,method=:nmf,
            delta=50ms,maxiter=100,tol=1e-3)

Ct = track(C,tc=1s)
rplot(Ct)

aCt = abs.(Ct);
thresh = mean(abs,Ct) + 1.5std(Ct);
tCt = copy(Ct);
tCt[aCt .< thresh] .= 0;
rplot(tCt)

tracks = any(abs.(Ct) .> thresh,(2,3))
df = DataFrame(y=Int.(vec(tracks)),
               x=repeat(ustrip(times(C)),outer=3),
               track=repeat(1:3,inner=ntimes(C)))
R"""
ggplot($df,aes(x,y + track*0.05,color=factor(track))) +
  scale_color_brewer(palette='Set1') +
  geom_point()
"""

function sources_active(window,thresh)
  min.(2,sum(any(window .> thresh,(1,2,3))))
end

function source_count_by_C(C;thresh_scale=1.2,window=0.25s,delta=0.1s,
                           buildup=1s)
  thresh = mean(C) + thresh_scale*std(C)
  map_window(C[buildup .. last(times(C))],window,delta) do win
    sources_active(win,thresh)
  end
end

counts = source_count_by_C(Ct,thresh_scale=2);
lplot(counts)

percept_lengths(counts)

