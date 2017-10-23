import Unitful: ms, Hz, kHz, s
using DataFrames
using Feather

include("model.jl")
include("stim.jl")

if !isdefined(:model)

end

fs = 8000
tone_len = 60ms

model = Model("/Users/davidlittle/Data")
a_stim = run_spect(model,tone(500Hz,tone_len,fs))
b_stim = run_spect(model,tone(500Hz * 2^(6/12),tone_len,fs))
stim = run_spect(model,aba(tone_len,10,fs,500Hz,6))

using Plots
# plotlyjs()
# gr()
order = unit_ordering_by(model.layer1,a_stim,b_stim)

function collect_responses(model,stim,taus=1:3)
  l1,l2,l3 = run(model,stim,return_all=true)

  indices = CartesianRange(size(l1))

  d = DataFrame(response = l1[:],
                time = [ii[1]*delta_t(model.layer1) for ii in indices],
                unit = [ii[1] for ii in indices],
                layer = "l1")

  for tau in taus
    di = DataFrame(response = l2[tau][:],
                   time = [ii[1]*delta_t(model.layer2) for ii in indices],
                   unit = [ii[1] for ii in indices],
                   layer = "l2_tau$(tau)")
    d = vcat(d,di)

    di = DataFrame(response = l3[tau][:],
                   time = [ii[1]*delta_t(model.layer2) for ii in indices],
                   unit = [ii[1] for ii in indices],
                   layer = "l3_tau$(tau)")
    d = vcat(d,di)
  end
  d
end

model = Model("/Users/davidlittle/Data",
              l1params = LayerParams(c_a=1,c_mi=1),
              l2params = LayerParmas(c_a=0,c_mi=0))
df = collect_responses(model,stim)
filename = "../../data/responses_l1_adapt.feather"
Feather.write(filename,df)
println("Wrote results to $filename.")

# model = Model("/Users/davidlittle/Data",l1params = L1Params(c_mi = 0))
# l1_A = run(model.layer1,stim)
# heatmap((indices(l1,1)-4) * (frame_length(model)/s),linspace(0,1,size(l1,2)),
#         l1[5:end,[order[1:20]; order[end-20:end]]]',
#         subplot = 3)

# plot!(model.spect,stim,subplot=4,c=:speed)

# filename = ("../../plots/adapt_visual_"*
#             Dates.format(now(),"yyyy-mm-dd")*".png")
# png(filename)

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
