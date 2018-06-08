# TODO: task - find a reasonable parameter set
# something seems to be different with the new moderating
# parametesr
using RCall
R"library(ggplot2)"
using ProgressMeter
include("count_lengths.jl")

dir = "../../../plots/run_$(Date(now()))"
isdir(dir) || mkdir(dir)

stim = ab(120ms,120ms,1,25,500Hz,6) |> normpower |> amplify(-10dB)
stim_resp = cortical(audiospect(stim,progressbar=false),progressbar=false,
                     scales=cycoct.*2.0.^linspace(-1,2,9))

methods = Dict(:threshold =>
               x->source_count_by_threshold(x,window=1s,delta=0.25s,
                                           cutoff=2cycoct,buildup=1s))
params = Dict(:c_σ=>0.3, :τ_σ=>500ms, :c_m=>30, :τ_m=>350ms, :W_m_σ => 10,
              :c_a=>10, :τ_a=>3s, :c_x=>3.0, :τ_x=>300ms, :c_n=>15, :τ_n=>10s)
# y,a,m = bistable_scales(stim_resp,params,intermediate_results=true)
# rplot(y[1s .. 10s])
# rplot(a[1s .. 10s])
# rplot(m[1s .. 10s])

# all_rows = []
@showprogress for i in 1:500
  rows = count_lengths_helper(stim_resp,methods,0,Date(now()),params)
  push!(all_rows,rows)
end

final_rows = vcat(all_rows...)
jldopen("single_histogram.jld2","w") do file
  file["rows"] = final_rows
end
alert("DONE!!!!")

df = DataFrame(length=map(x->x.length,final_rows),
               stimulus=map(x->x.stimulus,final_rows),
               error_code=map(x->x.error_code,final_rows))

R"""
library(cowplot)
ggplot(subset($df,error_code == 0),aes(x=length/mean(length))) +
  geom_histogram(bins=20) + ylab('') + xlab('') +
  coord_cartesian(xlim=c(-1,6))
ggsave($(joinpath(dir,"histogram.pdf")))
"""

