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

# TODO: what happens when we get ride of the scale filter

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect) #,scales=[2,4,8],rates=2.^(1:5))
tempc = TCAnalysis(cort,10)

ab0_1_async = @> ab(60ms,240ms,120ms,20,500Hz,0.1) attenuate(10)
ab1_async = @> ab(60ms,240ms,120ms,20,500Hz,1) attenuate(10)

ab6_sync = @> ab(60ms,240ms,0ms,20,500Hz,6) attenuate(10)
ab6_async = @>(ab(60ms,240ms,120ms,20,500Hz,6),attenuate(10))

ab12_async = @> ab(60ms,240ms,120ms,20,500Hz,12) attenuate(10)

ts, = tempc(ab0_1_async); ts[1] / ts[3]
ts, = tempc(ab1_async); ts[1] / ts[3]
ts, = tempc(ab6_async); ts[1] / ts[3]
ts, = tempc(ab12_async); ts[1] / ts[3]

ts, = tempc(ab6_sync); ts[1] / ts[3]

# TODO: show the rate only implementation here

cort = CorticalModel(spect,scales=[NaN],rates=2.^(1:5))
tempc = TCAnalysis(cort,10)

ts, = tempc(ab0_1_async); ts[1] / ts[3]
ts, = tempc(ab1_async); ts[1] / ts[3]
ts, = tempc(ab6_async); ts[1] / ts[3]
ts, = tempc(ab12_async); ts[1] / ts[3]

ts, = tempc(ab6_sync); ts[1] / ts[3]
