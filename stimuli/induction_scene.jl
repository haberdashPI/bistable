push!(LOAD_PATH,"../models/little/packages")
using AuditoryModel
using Sounds
include("abas.jl")

# context1 = aba(freq=2kHz,df=0,tone_length=75ms,spacing=50ms,repeat=32)
test = aba(freq=500Hz,df=6,tone_length=75ms,spacing=50ms,repeat=32)
context1 = aba_band(freq=2.5kHz,width=1,df=0,tone_length=75ms,spacing=50ms,repeat=32)
context2 = aba(freq=2kHz,df=missing,tone_length=75ms,spacing=50ms,repeat=32)

save("context1.wav",[amplify(context1,-10dB); amplify(test,-20)])
save("context2.wav",[amplify(context2,-10dB); amplify(test,-20)])

baseline = aba(freq=500Hz,df=12,tone_length=75ms,spacing=50ms,repeat=64)
broad = aba_band(freq=500Hz,width=0.5,df=12,tone_length=75ms,spacing=50ms,repeat=64)
pullsed = aba_pulse(freq=500Hz,pulse=2,pulse_spacing=5ms,df=6,tone_length=75ms,spacing=50ms,repeat=64)

save("baseline.wav",amplify(baseline,-10dB))
save("broad.wav",amplify(broad,-10dB))
save("pullsed.wav",amplify(pullsed,-10dB))

rplot(audiospect(pullsed))

# alter_4 = [aba(freq=500Hz,df=12,tone_length=75ms,spacing=50ms,repeat=4);
#            aba(freq=500Hz,df=1,tone_length=75ms,spacing=50ms,repeat=4)]

# alter_8 = [aba(freq=500Hz,df=12,tone_length=75ms,spacing=50ms,repeat=8);
#            aba(freq=500Hz,df=1,tone_length=75ms,spacing=50ms,repeat=8)]

# save("aba_alter_4.wav",vcat(Iterators.repeated(alter_4,4)...))
# save("aba_alter_8.wav",vcat(Iterators.repeated(alter_8,2)...))
