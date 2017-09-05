import Unitful: ms, Hz, kHz, s
using Plots; plotlyjs()

include("model.jl")
include("stim.jl")
include("audio_spect.jl")

audio_spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5")

if !isdefined(:model)
  model = Model("/Users/davidlittle/Data/model.h5")
end

fs = 8000
tone_len = 60ms
a_stim = run(audio_spect,tone(500Hz,tone_len,fs))
b_stim = run(audio_spect,tone(500Hz * 2^(6/12),tone_len,fs))

stim = run(audio_spect,aba(tone_len,10,fs,500Hz,6))

l1 = run(model.layer1,Float32.(stim))
a_match = resemblance(model.layer1,l1,Float32.(a_stim[3:6,:]))
b_match = resemblance(model.layer1,l1,Float32.(b_stim[3:6,:]))

# gr()
plot(indices(a_match,1) * 1/fs * frame_length(audio_spect) * frame_length(model),
     [a_match b_match],label = ["A" "B"],
     subplot = 1,layout = grid(2,1,heights = [0.8,0.2]))

plot!(audio_spect,stim,subplot=2,c=:speed)

# filename = ("../../plots/ab_visual_"*
#               Dates.format(now(),"yyyy-mm-dd")*".png")
# savefig(filename)

# TODO: start using plausible stimulus: these are way too fast.
# TODO: test the plausible stimulus on the build-up and sycnhronization
# tasks.
# TODO: create the adaptation version of the model, possibly
# repeating frames to get a good "resolution".
