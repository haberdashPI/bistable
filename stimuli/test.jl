using TimedSound
# include("../models/little/units.jl")
include("warble.jl")

A1_ = @> warble(0.5kHz           ,0.5,0.8,3,4,20s) sound attenuate(10) ramp(50ms)
B1_ = @> warble(0.5kHz * 2^(1/12),0.5,0.8,3,4,20s) sound attenuate(10) ramp(50ms)
x = @> mix(A1_,B1_) attenuate(15)

@> mix(A1_,B1_) attenuate(15) play

B1__ = @> fadeto(silence(75ms + 10ms),B1_,10ms)
@> mix(A1_,B1__) attenuate(15) play


A4_ = @> warble(0.5kHz           ,0.5,0.8,3,4,20s) sound attenuate(10) ramp(50ms)
B4_ = @> warble(0.5kHz * 2^(3/12),0.5,0.8,3,4,20s) sound attenuate(10) ramp(50ms)
xs = @> mix(A4_,B4_) attenuate(15)

B4__ = @> fadeto(silence(125ms + 10ms),B4_,10ms)
xa = @> mix(A4_[75ms .. ends],B4__[125ms .. ends]) attenuate(15)

A4off_ = @> warble(0.5kHz           ,0.5,0.8,3,5,20s) sound attenuate(10) ramp(50ms)
xo = @> mix(A4off_,B4__) attenuate(15)

A4off__ = @> warble(0.5kHz           ,0.5,0.8,3,7,20s) sound attenuate(10) ramp(50ms)
xoo = @> mix(A4off__,B4__) attenuate(15)

Ad_ = @> warble(2kHz           ,0.5,0.8,1,10,20s) sound attenuate(10) ramp(50ms)
xd = @> mix(Ad_,B4__) attenuate(15)


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
width = 1/12
a = @> noise(100ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(50ms)
freq = 1.5kHz
width = 6/12
b = @> noise(100ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(50ms)

x = vcat(Iterators.repeated([a; b],100)...)


freq = 1kHz
width = 8/12
a = @> noise(25ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(75ms)]
freq = 1.5kHz
width = 2/12
b = @> noise(25ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(75ms)]

x = vcat(Iterators.repeated([a; b],100)...)



freq = 1kHz
width = 8/12
a = @> noise(25ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(75ms)]
freq = 1.5kHz
width = 2/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(10ms)]

x = vcat(Iterators.repeated([a; b],100)...)



freq = 1kHz
width = 12/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = 1.5kHz
width = 1/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x = vcat(Iterators.repeated([a; b],100)...)


freq = 1kHz
width = 12/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = freq*2^(0/12)
width = 1/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x11 = vcat(Iterators.repeated([a; b],100)...)

freq = 1kHz
width = 8/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = freq*2^(0/12)
width = 4/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x4 = vcat(Iterators.repeated([a; b],100)...)


freq = 1kHz
width = 4/12
a = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
a = [a; silence(30ms)]
freq = freq*2^(0/12)
width = 4/12
b = @> noise(90ms) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2)) ramp(5ms)
b = [b; silence(30ms)]

x0 = vcat(Iterators.repeated([a; a],100)...)

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
b = wobble(2^(3/12)*1kHz,depth,4Hz,π,20s)

@> mix(a,b) attenuate(15) play


depth = 0.5
a = wobble(1kHz,depth,4Hz,0,20s)
b = wobble(2^(3/12)*1kHz,depth,4Hz,π,20s)

@> mix(a,b) attenuate(15) play



depth = 1
a = wobble(1kHz,depth,4Hz,0,20s)
b = wobble(2^(3/12)*1kHz,depth,4Hz,π,20s)

@> mix(a,b) attenuate(15) play

depth = 25
a = wobble(1kHz,depth,4Hz,0,20s)
b = wobble(2^(3/12)*1kHz,depth,4Hz,π,20s)

@> mix(a,b) attenuate(15) play

# seems more like the above adjusts the point of bistability rather
# than eliminating it

# what do these variations in dimensions buy us in understanding? what can we
# learn from it about the model?

# how does this work when we vary along multiple dimensions and approach
# a point of bistability? can we learn something from that?
