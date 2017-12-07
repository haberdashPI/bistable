# NOTES from meeting
#=
The key point for getting the coherence working is to make sure
that there are time delays that allow stimuli across time to be
compared.

tasks:

1. figure out if the support vectors themselves are informative
enough of the number of objects, (probably use a different set of rates),
why does the A stimulus disappear? (is there something wrong with the
cortical model representation?)

2. look through papers on envelope detection, think through a possible
eeg experiment using this for merve's paper, and/or think through
potential behavioral studies
=#

include("units.jl")
include("tempc.jl")
include("stim.jl")

# TODO: try the current approach but
# limit to a single temporal rate
# to try and reporduce the results from
# the 2009 paper

# (may require looking at a rate only filter)

setup_sound(sample_rate=8kHz)


# TODO: create plots for the below
# results

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)

ab0_1_async = @> ab(60ms,60ms,1,20,500Hz,0.1) attenuate(10)
ab1_async = @> ab(60ms,60ms,1,20,500Hz,1) attenuate(10)

ab6_sync = @> ab(60ms,60ms,0,20,500Hz,6) attenuate(10)
ab6_async = @>(ab(60ms,60ms,1,20,500Hz,6),attenuate(10))

ab12_async = @> ab(60ms,60ms,1,20,500Hz,12) attenuate(10)
ab12_async_240 = @> ab(4*60ms,4*60ms,1,20,500Hz,12) attenuate(10)

tempc = TCAnalysis(CorticalModel(spect),10,4s)
C = tempc(ab12_async_240)
rplot(tempc,C,oddonly=true)

dir = "../../plots/meeting_2017_12_07"
p = Array{Any}(16)
########################################
# plot 1 - compare models on an unambiguous stimulus

# TODO: clean up the graph

tempc = TCAnalysis(CorticalModel(spect),10,4s)
C = tempc(ab12_async_240)
p[1] = rplot(tempc,C[end-10].λ)

tempc = TCAnalysis(CorticalModel(spect,scales=[NaN],rates=2.^(1:5)),10,4s)
C = tempc(ab12_async_240)
p[2] = rplot(tempc,C[end-10].λ)

tempc = TCAnalysis(CorticalModel(spect,scales=[NaN],
                                 rates=[-2.0.^(1:5);2.^(1:5)]),10,4s)
C = tempc(ab12_async_240)
p[3] = rplot(tempc,C[end-10].λ)

tempc = TCAnalysis(CorticalModel(spect,rates=2.^(1:5)),10,4s)
C = tempc(ab12_async_240)
p[4] = rplot(tempc,C[end-10].λ)

R"""
library(cowplot)

p = plot_grid($(p[1]) + ggtitle('full model'),
  $(p[4]) + ggtitle('+rates only'),
  $(p[2]) + ggtitle('no scales, +rates only'),
  $(p[3]) + ggtitle('no scales, all rates'),
  align="h",ncol=2)

save_plot(paste($dir,"/eigeval_model_varations.pdf",sep=''),
  p,ncol=2,nrow=2,base_aspect_ratio=1)
"""

########################################
# plot 2 - plot of unambiguous stimulus components
p = Array{Any}(2)

tempc = TCAnalysis(CorticalModel(spect),10,4s)
C = tempc(ab12_async_240)
p[1] = rplot(tempc,C[end-10],n=6)

tempc = TCAnalysis(CorticalModel(spect,scales=[NaN],rates=2.^(1:5)),10,4s)
C = tempc(ab12_async_240)
p[2] = rplot(tempc,C[end-10],n=6)

R"""
library(cowplot)

p = plot_grid(
  $(p[1]) + ggtitle('Components, Full Model'),
  $(p[2]) + ggtitle('Components, no scales, +rates only'),
  align="h")

save_plot(paste($dir,"/components_model_variations.pdf",sep=''),
  p,ncol=2,nrow=1,base_aspect_ratio=1.3,base_width=8)
"""
########################################
# plot 3 - eigen value plots of varying stimulus conditions for full model

p = Array{Any}(6)

tempc = TCAnalysis(CorticalModel(spect),10,4s)
C = tempc(ab12_async_240)
p[1] = rplot(tempc,C[end-10].λ)

C = tempc(ab0_1_async)
p[2] = rplot(tempc,C[end-10].λ)

C = tempc(ab1_async)
p[3] = rplot(tempc,C[end-10].λ)

C = tempc(ab6_async)
p[4] = rplot(tempc,C[end-10].λ)

C = tempc(ab6_sync)
p[5] = rplot(tempc,C[end-10].λ)


C = tempc(ab12_async)
p[6] = rplot(tempc,C[end-10].λ)

R"""
library(cowplot)

p = plot_grid(
  $(p[1]) + ggtitle('12 st, 240ms unit'),
  $(p[2]) + ggtitle('0.1 st, 60ms unit'),
  $(p[3]) + ggtitle('1 st, 60ms unit'),
  $(p[4]) + ggtitle('6 st, 60ms unit'),
  $(p[5]) + ggtitle('6 st synchronous, 60ms unit'),
  $(p[6]) + ggtitle('12 st, 60ms unit'),
  align="h",nrow=2)

save_plot(paste($dir,"/eigenval_stimulus_variations.pdf",sep=''),
  p,ncol=3,nrow=2,base_aspect_ratio=1.3,base_width=6)
"""

########################################
# plot 3 - plot of eigenvalue time course

tempc = TCAnalysis(CorticalModel(spect),10,4s)
p[1] = rplot(tempc,ab12_async_240,n=6,oddonly=true)
p[2] = rplot(tempc,ab12_async,n=6,oddonly=true)
p[3] = rplot(tempc,ab6_sync,n=6,oddonly=true)

R"""
library(cowplot)

p = plot_grid(
  $(p[1]) + ggtitle('12 st, 240ms unit'),
  $(p[2]) + ggtitle('12 st, 60ms unit'),
  $(p[3]) + ggtitle('6 st synchronous, 60ms unit'),
  align="h",ncol=1)

save_plot(paste($dir,"/timecourse.pdf",sep=''),p,
          ncol=1,nrow=3,base_aspect_ratio=2.4,base_width=9)
"""

########################################
# plot 4 - plot of eigencomponents of different stimuli

tempc = TCAnalysis(CorticalModel(spect),10,4s)
C = tempc(ab12_async_240)
p[1] = rplot(tempc,C[end-10],n=6)

C = tempc(ab12_async)
p[2] = rplot(tempc,C[end-10],n=6)

C = tempc(ab6_sync)
p[3] = rplot(tempc,C[end-10],n=6)


R"""
library(cowplot)

p = plot_grid(
  $(p[1]) + ggtitle('12 st, 240ms unit'),
  $(p[2]) + ggtitle('12 st, 60ms unit'),
  $(p[3]) + ggtitle('6 st synchronous, 60ms unit'),
  align="h",ncol=1)

save_plot(paste($dir,"/components_stimulus_variations.pdf",sep=''),p,
          ncol=1,nrow=3,base_aspect_ratio=1.3,base_width=6)
"""

##################################################
# plot 5 - of eigen values calculated seperately for each rate
tempc = TCAnalysis(CorticalModel(spect),10,4s,split_rates=true)
C12 = tempc(ab12_async)
C6 = tempc(ab6_async)
C1 = tempc(ab1_async)
C01 = tempc(ab0_1_async)

N = length(rates(tempc.upstream))
df = vcat(DataFrame(level = vcat((C12[i][end-10].λ[[1,3]] for i in 1:N)...),
                    st = string(12)*"st",
                    rate = repeat(rates(tempc.upstream),inner=2),
                    component = repeat([1,3],outer=N)),

          DataFrame(level = vcat((C6[i][end-10].λ[[1,3]] for i in 1:N)...),
                    st = "0"*string(6)*"st",
                    rate = repeat(rates(tempc.upstream),inner=2),
                    component = repeat([1,3],outer=N)),

          DataFrame(level = vcat((C1[i][end-10].λ[[1,3]] for i in 1:N)...),
                    st = "0"*string(1)*"st",
                    rate = repeat(rates(tempc.upstream),inner=2),
                    component = repeat([1,3],outer=N)),

          DataFrame(level = vcat((C01[i][end-10].λ[[1,3]] for i in 1:N)...),
                    st = string(0.1)*"st",
                    rate = repeat(rates(tempc.upstream),inner=2),
                    component = repeat([1,3],outer=N)))


R"""
library(ggplot2)
library(cowplot)

p = ggplot($df,aes(x=factor(round(rate)),y=level,fill=factor(component))) +
  geom_bar(stat='identity',position='dodge') + facet_grid(st~.) +
  scale_fill_brewer(palette='Set1',name='Component') +
  ggtitle('Eigenvalues computed separately for each rate')

save_plot(paste($dir,"/eigvals_by_rate.pdf",sep=''),
  p,base_aspect_ratio=2.4,base_width=8) + xlab('rate (Hz)')
"""


# ii = CartesianRange(size(C.λ))
# slicei(i) = map(ii -> ii[i],ii)
# df = DataFrame(level = vec(C.λ),
#                time = vec(ustrip(slicei(1) * Δt(spect))),
#                component = vec(slicei(2)))

# R"""
# library(ggplot2)

# ggplot(subset($df,(component %% 2) == 0),
#   aes(x=time,y=level,color=factor(component),group=component)) +
#   geom_line()
# """

# λ, = tempc(ab6_sync); λ[1] / λ[3]

# λ,ϕ,v = tempc(ab12_async)
# rplot(spect,Diagonal(λ)*squeeze(mean(ϕ,2),2))
# rplot(spect,squeeze(mean(v,2),2))
# rplot(CorticalModel(spect,scales=[NaN],rates=rates(cort)),
#       reshape(v,size(v,1:2...)...,1,:))
# λ[1] / λ[3]

# λ,ϕ,v = tempc(ab1_async)
# rplot(spect,Diagonal(λ)*squeeze(mean(ϕ,2),2))
# rplot(spect,squeeze(mean(v,2),2))
# rplot(CorticalModel(spect,scales=[NaN],rates=rates(cort)),
#       reshape(v,size(v,1:2...)...,1,:))
# λ[1] / λ[3]

# cort = CorticalModel(spect,scales=[NaN],rates=2.^(1:5))
# tempc = TCAnalysis(cort,10,method=:batch)

# λ, = tempc(ab0_1_async); λ[1] / λ[2]
# λ, = tempc(ab1_async); λ[1] / λ[2]
# λ, = tempc(ab6_async); λ[1] / λ[2]
# λ, = tempc(ab12_async); λ[1] / λ[2]

# λ, = tempc(ab6_sync); λ[1] / λ[2]

# λ,ϕ,v = tempc(ab12_async)
# λ,ϕ,v = tempc(ab1_async)

# # can we can add doubling by adding negative rates? nope
# cort = CorticalModel(spect,scales=[NaN])
# tempc = TCAnalysis(cort,10,method=:batch)

# λ, = tempc(ab0_1_async); λ[1] / λ[2]
# λ, = tempc(ab1_async); λ[1] / λ[2]
# λ, = tempc(ab6_async); λ[1] / λ[2]
# λ, = tempc(ab12_async); λ[1] / λ[2]

# λ, = tempc(ab6_sync); λ[1] / λ[2]

# # can we can remove doubling by removing negative rates? nope
# cort = CorticalModel(spect,rates=2.^(1:5))
# tempc = TCAnalysis(cort,10,method=:batch)
# λ, = tempc(ab0_1_async)
# λ, = tempc(ab12_async)


# cort = CorticalModel(spect) #,scales=[2,4,8],rates=2.^(1:5))
# tempc = TCAnalysis(cort,10,method=:ipca)
# λ, = tempc(ab0_1_async); λ[end,end] / λ[end,end-3]
# λ, = tempc(ab1_async);  λ[end,end] / λ[end,end-3]
# λ, = tempc(ab6_async);  λ[end,end] / λ[end,end-3]
# λ, = tempc(ab12_async); λ[end,end] / λ[end,end-3]

# λ, = tempc(ab6_sync);  λ[end,end-1] / λ[end,end-3]

# λ,ϕ,v = tempc(ab12_async)
# rplot(spect,Diagonal(λ[end,:])*squeeze(mean(ϕ[end,:,:,:],2),2))
# rplot(spect,v)
# rplot(CorticalModel(spect,scales=scales(cort),rates=1:size(ϕ,2)),
#       ϕ[:,:,:,:] .* λ[1:1,:,1:1,1:1])
# λ[end,end] / λ[end,end-3]

# df = DataFrame(time=times(spect,λ),
#                ratio=λ[:,end] ./ λ[:,end-3])
# R"""
# ggplot($df,aes(x=time,y=ratio)) + geom_line() + ggtitle("Delta 12 s.t.")
# """

# λ,ϕ,v = tempc(ab1_async)
# rplot(spect,Diagonal(λ[end,:])*squeeze(mean(ϕ[end,:,:,:],2),2))
# rplot(spect,v)
# rplot(CorticalModel(spect,scales=scales(cort),rates=1:size(ϕ,2)),
#       ϕ[:,:,:,:] .* λ[1:1,:,1:1,1:1])
# λ[end,end] / λ[end,end-3]

# df = DataFrame(time=times(spect,λ),
#                ratio=λ[:,end] ./ λ[:,end-3])
# R"""
# ggplot($df,aes(x=time,y=ratio)) + geom_line() + ggtitle("Delta 12 s.t.")
# """

# λ1,ϕ,v = tempc(ab1_async)
# λ12,ϕ,v = tempc(ab12_async)

# df = DataFrame(time=[times(spect,λ); times(spect,λ)],
#                ratio=[λ1[:,end] ./ λ1[:,end-3];
#                       λ12[:,end] ./ λ12[:,end-3]],
#                st = [fill(1,size(λ1,1)); fill(12,size(λ1,1))])
# R"""
# ggplot($df,aes(x=time,y=ratio,color=factor(st))) + geom_line()
# """


# ab1_asyncL = @> ab(120ms,480ms,240ms,20,500Hz,1) attenuate(10)
# ab12_asyncL = @> ab(120ms,480ms,240ms,20,500Hz,12) attenuate(10)
# λ1,ϕ,v = tempc(ab1_asyncL)
# λ12,ϕ,v = tempc(ab12_asyncL)

# df = DataFrame(time=[times(spect,λ1); times(spect,λ1)],
#                ratio=[λ1[:,end] ./ λ1[:,end-3];
#                       λ12[:,end] ./ λ12[:,end-3]],
#                st = [fill(1,size(λ1,1)); fill(12,size(λ1,1))])
# R"""
# ggplot($df,aes(x=time,y=ratio,color=factor(st))) + geom_line()
# """

# TODO: improve on visualizaiton of component amplitude -
# TODO: improve on visualization of components themselves,
# these plots remain relatively unclear.

# TODO: start implementing variations on
# adaptation-MI in the three levels
# for this I want the strengths of the
# eigen values to adapt and mutually inhibit
# goal for thursday: show that MI of lower layers alone
# do not suffice
