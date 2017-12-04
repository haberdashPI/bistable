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

ab0_1_async = @> ab(60ms,240ms,120ms,20,500Hz,0.1) attenuate(10)
ab1_async = @> ab(60ms,240ms,120ms,20,500Hz,1) attenuate(10)

ab6_sync = @> ab(60ms,240ms,0ms,20,500Hz,6) attenuate(10)
ab6_async = @>(ab(60ms,240ms,120ms,20,500Hz,6),attenuate(10))

ab12_async = @> ab(60ms,240ms,120ms,20,500Hz,12) attenuate(10)

cort = CorticalModel(spect) #,scales=[2,4,8],rates=2.^(1:5))
tempc = TCAnalysis(cort,10,method=:batch)
λ, = tempc(ab0_1_async); λ[1] / λ[3]
λ, = tempc(ab1_async); λ[1] / λ[3]
λ, = tempc(ab6_async); λ[1] / λ[3]
λ, = tempc(ab12_async); λ[1] / λ[3]

λ, = tempc(ab6_sync); λ[1] / λ[3]

λ,ϕ,v = tempc(ab12_async)
rplot(spect,Diagonal(λ)*squeeze(mean(ϕ,2),2))
rplot(spect,squeeze(mean(v,2),2))
rplot(CorticalModel(spect,scales=[NaN],rates=rates(cort)),
      reshape(v,size(v,1:2...)...,1,:))
λ[1] / λ[3]

λ,ϕ,v = tempc(ab1_async)
rplot(spect,Diagonal(λ)*squeeze(mean(ϕ,2),2))
rplot(spect,squeeze(mean(v,2),2))
rplot(CorticalModel(spect,scales=[NaN],rates=rates(cort)),
      reshape(v,size(v,1:2...)...,1,:))
λ[1] / λ[3]

cort = CorticalModel(spect,scales=[NaN],rates=2.^(1:5))
tempc = TCAnalysis(cort,10,method=:batch)

λ, = tempc(ab0_1_async); λ[1] / λ[2]
λ, = tempc(ab1_async); λ[1] / λ[2]
λ, = tempc(ab6_async); λ[1] / λ[2]
λ, = tempc(ab12_async); λ[1] / λ[2]

λ, = tempc(ab6_sync); λ[1] / λ[2]

λ,ϕ,v = tempc(ab12_async)
λ,ϕ,v = tempc(ab1_async)

# can we can add doubling by adding negative rates? nope
cort = CorticalModel(spect,scales=[NaN])
tempc = TCAnalysis(cort,10,method=:batch)

λ, = tempc(ab0_1_async); λ[1] / λ[2]
λ, = tempc(ab1_async); λ[1] / λ[2]
λ, = tempc(ab6_async); λ[1] / λ[2]
λ, = tempc(ab12_async); λ[1] / λ[2]

λ, = tempc(ab6_sync); λ[1] / λ[2]

# can we can remove doubling by removing negative rates? nope
cort = CorticalModel(spect,rates=2.^(1:5))
tempc = TCAnalysis(cort,10,method=:batch)
λ, = tempc(ab0_1_async)
λ, = tempc(ab12_async)


cort = CorticalModel(spect) #,scales=[2,4,8],rates=2.^(1:5))
tempc = TCAnalysis(cort,10,method=:ipca)
λ, = tempc(ab0_1_async); λ[end,end] / λ[end,end-3]
λ, = tempc(ab1_async);  λ[end,end] / λ[end,end-3]
λ, = tempc(ab6_async);  λ[end,end] / λ[end,end-3]
λ, = tempc(ab12_async); λ[end,end] / λ[end,end-3]

λ, = tempc(ab6_sync);  λ[end,end-1] / λ[end,end-3]

λ,ϕ,v = tempc(ab12_async)
rplot(spect,Diagonal(λ[end,:])*squeeze(mean(ϕ[end,:,:,:],2),2))
rplot(spect,v)
rplot(CorticalModel(spect,scales=scales(cort),rates=1:size(ϕ,2)),
      ϕ[:,:,:,:] .* λ[1:1,:,1:1,1:1])
λ[end,end] / λ[end,end-3]

df = DataFrame(time=times(spect,λ),
               ratio=λ[:,end] ./ λ[:,end-3])
R"""
ggplot($df,aes(x=time,y=ratio)) + geom_line() + ggtitle("Delta 12 s.t.")
"""

λ,ϕ,v = tempc(ab1_async)
rplot(spect,Diagonal(λ[end,:])*squeeze(mean(ϕ[end,:,:,:],2),2))
rplot(spect,v)
rplot(CorticalModel(spect,scales=scales(cort),rates=1:size(ϕ,2)),
      ϕ[:,:,:,:] .* λ[1:1,:,1:1,1:1])
λ[end,end] / λ[end,end-3]

df = DataFrame(time=times(spect,λ),
               ratio=λ[:,end] ./ λ[:,end-3])
R"""
ggplot($df,aes(x=time,y=ratio)) + geom_line() + ggtitle("Delta 12 s.t.")
"""

λ1,ϕ,v = tempc(ab1_async)
λ12,ϕ,v = tempc(ab12_async)

df = DataFrame(time=[times(spect,λ); times(spect,λ)],
               ratio=[λ1[:,end] ./ λ1[:,end-3];
                      λ12[:,end] ./ λ12[:,end-3]],
               st = [fill(1,size(λ1,1)); fill(12,size(λ1,1))])
R"""
ggplot($df,aes(x=time,y=ratio,color=factor(st))) + geom_line()
"""


ab1_asyncL = @> ab(120ms,480ms,240ms,20,500Hz,1) attenuate(10)
ab12_asyncL = @> ab(120ms,480ms,240ms,20,500Hz,12) attenuate(10)
λ1,ϕ,v = tempc(ab1_asyncL)
λ12,ϕ,v = tempc(ab12_asyncL)

df = DataFrame(time=[times(spect,λ1); times(spect,λ1)],
               ratio=[λ1[:,end] ./ λ1[:,end-3];
                      λ12[:,end] ./ λ12[:,end-3]],
               st = [fill(1,size(λ1,1)); fill(12,size(λ1,1))])
R"""
ggplot($df,aes(x=time,y=ratio,color=factor(st))) + geom_line()
"""

# TODO: improve on visualizaiton of component amplitude -
# TODO: improve on visualization of components themselves,
# these plots remain relatively unclear.

# TODO: start implementing variations on
# adaptation-MI in the three levels
# for this I want the strengths of the
# eigen values to adapt and mutually inhibit
# goal for thursday: show that MI of lower layers alone
# do not suffice
