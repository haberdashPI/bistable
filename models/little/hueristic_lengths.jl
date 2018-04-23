push!(LOAD_PATH,"packages")
using AuditoryModel
using ProgressMeter
using DataFrames
using AuditoryCoherence
using AxisArrays
using RCall
include("util/stim.jl")
include("util/lengths.jl")

R"library(ggplot2)"
# R"library(cowplot)"
# quartz() = R"quartz()"
dir = "../../plots/run_2018_04_02"
isdir(dir) || mkdir(dir)

window(fn,win,xs::AbstractVector) = [fn(xs[clamp.(i:i+win,0,end)])
                                     for i in 1:win:length(xs)]

function percept_lengths(csa,windowlen,buildup)
  threshold = mean(abs,csa)

  two_stream_scale = 2cycoct

  onestream = sum(max.(0,abs.(csa[:,0cycoct .. two_stream_scale,:]) .-
                       threshold),
                  [axisdim(csa,Axis{:scale}) axisdim(csa,Axis{:freq})])
  twostream = sum(max.(0,abs.(csa[:,two_stream_scale .. 4cycoct,:]) .-
                       threshold),
                  [axisdim(csa,Axis{:scale}) axisdim(csa,Axis{:freq})])

  onemean = window(mean,windowlen,vec(onestream))
  twomean = window(mean,windowlen,vec(twostream))
  N = length(onemean)
  ts = times(csa)[1:windowlen:ntimes(csa)]

  onemean ./= maximum(onemean[ts .> buildup])
  twomean ./= maximum(twomean[ts .> buildup])

  len,stim = findlengths(onemean .> twomean)
  len .* (ts[2] - ts[1]),stim
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
  len,stim = percept_lengths(y,20,0.5s)
  push!(lengths,len)
  push!(stimuli,stim)
end
df = DataFrame(length = ustrip.(vcat(lengths...)), stimuli = vcat(stimuli...))

R"""
ggplot($df,aes(x=length,fill=stimuli)) +
  geom_histogram(bins=30) + xlab('Percept Length (s)') +
  scale_fill_brewer(palette='Set1',name='Streaming',
                     labels=c('One Stream','Two Streams')) +
  theme_classic()
ggsave($(joinpath(dir,"histogram_bythreshold.pdf")))
"""
