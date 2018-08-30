push!(LOAD_PATH,"packages")
using AuditoryModel
using ProgressMeter
using DataFrames
using AuditoryCoherence
using AxisArrays
using RCall

include("util/stim.jl")
include("util/peaks.jl")
include("util/lengths.jl")
include("util/biapply.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_04_12"
isdir(dir) || mkdir(dir)

function percept_lengths(csa,minlen=1s)
  counts = source_count_by_peaks(csa)
  thresh = ustrip(uconvert(s,minlen/Δt(counts)))
  lens,vals = findlengths(counts)
  slens = lens * ustrip(uconvert(s,Δt(counts)))

  mergelengths(slens,vals,ustrip(uconvert(s,minlen)))
end


noise_params = Dict(:τ_σ => 500ms,:c_σ => 0.3)
adapt_params = Dict(:c_m=>5,:τ_m=>350ms,:c_e=>1.5,:τ_e => 300ms,:c_a=>30,
                    :τ_a=>3s,:α=>3.0)
cortical_params = Dict(:scales=>cycoct.*round.(2.0.^linspace(-1,2,18),1))

x = ab(120ms,120ms,1,50,500Hz,6) |> normpower |> amplify(-10dB)

lengths = []
stimuli = []
@showprogress for i = 1:100
  y = bistable_scales(x,noise_params,adapt_params,cortical_params)
  len,stim = percept_lengths(y)
  push!(lengths,len)
  push!(stimuli,stim)
end
df = DataFrame(length = ustrip.(vcat(lengths...)), stimuli = vcat(stimuli...))

R"""
ggplot($df,aes(x=length,fill=factor(stimuli))) +
  geom_histogram(bins=30) + xlab('Percept Length (s)') +
  scale_fill_brewer(palette='Set1',name='Streaming',
                     labels=c('One Stream','Two Streams')) +
  theme_classic()
ggsave($(joinpath(dir,"histogram_bypeaks.pdf")))
"""
