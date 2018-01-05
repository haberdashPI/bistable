using TimedSound
include("../models/little/units.jl")
include("warble.jl")

# thoughts: these objects have no shapr boundaries
# which may lead to more ambguity

# problem: even the single sounds have some ambiguity to them
# I need something that clearly sounds like an object

# hypothesis - attention exagerates adapt/mi along
#   salient/and or top-down directed features
#   what is actually bistable in these
#   paradigms is the locus of attention
#   there are potentially lower-level loci but
#   these are transitory, attenion drives the ongoing
#   oscillation

# proposal: the oscilation occurs when there are equally salient/directed
# regions... but that doesn't entirely explain forced bistability,


# note, bistability occurs even for speech sounds, if they are repeated

# in general the fluctuations occur for repeated stimuli, where adaptation
# can take effect...


# what properties does a sound need to be a potentially bistable object-grouping?
# - it needs to have a different number of modes across rate/scale
# - these modes need to be relatively balanced
# - there need to be salient onsets associated with
#   the different numbered modes
#
# prediction, if the low level features play a role in bistability: a number of
# stimulus parameters would prevent bistability i.e. including less energy in
# one of the modes either by reducing the length or bandwidth of B
# would limit the ability of the A and B to separate

freq = 1kHz
width = 12/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = freq*2^(0/12)
width = 1/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x11 = vcat(Iterators.repeated([a; b],100)...)
save("approach1_split.wav",x11)

freq = 1kHz
width = 8/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = freq*2^(0/12)
width = 4/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x4 = vcat(Iterators.repeated([a; b],100)...)
save("approach1_bistable.wav",x4)

freq = 1kHz
width = 4/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = freq*2^(0/12)
width = 4/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x0 = vcat(Iterators.repeated([a; a],100)...)
save("approach1_fuse.wav",x0)

# prediction: more gradual modulations would reduce the bistability of the
# effect, but not alter the ambiguity

sig(x) = 0.5+0.5tanh(x)
wobble(freq,depth,dfreq,dphase,len) = audible(len) do t
  f = ustrip(TimedSound.inHz(freq))
  carrier = @. sin(2π*f*t)

  df = ustrip(TimedSound.inHz(dfreq))
  amp = @. sig(depth*sin(2π*df*t + dphase))
  carrier.*amp
end

depth = 0.2
a = wobble(1kHz,depth,4Hz,0,20s)
b = wobble(2^(1/12)*1kHz,depth,4Hz,π,20s)

x = @> mix(a,b) attenuate(15)
save("approach2_split.wav",x)

depth = 1
a = wobble(1kHz,depth,4Hz,0,20s)
b = wobble(2^(2/12)*1kHz,depth,4Hz,π,20s)

x = @> mix(a,b) attenuate(15)
save("approach2_bistable.wav",x)


depth = 25
a = wobble(1kHz,depth,4Hz,0,20s)
b = wobble(2^(2/12)*1kHz,depth,4Hz,π,20s)

x = @> mix(a,b) attenuate(15)
save("approach2_fuse.wav",x)

# seems more like the above adjusts the point of bistability rather
# than eliminating it

# though maybe, with more listening, it seems like it might also
# alter the sense I have of either: that is, I'm not
# as certain about my responses.

# what do these variations in dimensions buy us in understanding? what can we
# learn from it about the model?

# how does this work when we vary along multiple dimensions and approach
# a point of bistability? can we learn something from that?
