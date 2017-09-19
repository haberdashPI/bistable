import Unitful: ms, Hz, kHz, s

include("model.jl")
include("stim.jl")

if !isdefined(:model)
  model = Model("/Users/davidlittle/Data")
end

fs = 8000
tone_len = 60ms

a_stim = run_spect(model,tone(500Hz,tone_len,fs))
b_stim = run_spect(model,tone(500Hz * 2^(6/12),tone_len,fs))
stim = run_spect(model,aba(tone_len,10,fs,500Hz,6))
l1,l1_A = run(model.layer1,stim,return_adaptation=true)

using Plots
plotlyjs()
 
order = unit_ordering_by(model.layer1,a_stim,b_stim)

# TODO: trying to visualize the more "A-ish" and more "B-ish" responses
heatmap((indices(l1,1)-4) * frame_length(model),linspace(0,1,size(l1,2)),
        l1[5:end,[order[1:20]; order[end-20:end]]]',
        subplot = 1,layout = grid(3,1,heights = [0.4,0.4,0.2]))

heatmap!((indices(l1_A,1)-4) * (frame_length(model)/model.layer1.N_a),
         linspace(0,1,size(l1,2)),
         l1_A[5*model.layer1.N_a:end,[order[1:20]; order[end-20:end]]]',
         subplot = 2)

plot!(model.spect,stim,subplot=3,c=:speed)

# plot(indices(a_match,1) * 1/fs * frame_length(audio_spect) * frame_length(model),
#      [a_match b_match],label = ["A" "B"],
#      subplot = 1,layout = grid(2,1,heights = [0.8,0.2]))


# filename = ("../../plots/ab_visual_"*
#               Dates.format(now(),"yyyy-mm-dd")*".png")
# savefig(filename)

# TODO: start using plausible stimulus: these are way too fast.
# TODO: test the plausible stimulus on the build-up and sycnhronization
# tasks.
# TODO: create the adaptation version of the model, possibly
# repeating frames to get a good "resolution".
