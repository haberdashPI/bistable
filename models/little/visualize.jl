import Unitful: ms, Hz, kHz, s
using Plots; plotlyjs()
using Distributions
include("model.jl")
include("stim.jl")

model = Model("/Users/davidlittle/Data",
              l1params = LayerParams(c_mi = 0,c_a = 0,use_sig=true),
              l2params = LayerParams(c_mi = 0,c_a = 0,use_sig=false))

# NOTE: looking at output of single units makes basically zero sense. The
# response is dominated by the bias, meaning that units only effectively alter
# the output when multiple units are active.
# heatmap(sig.(model.layer1.w .+ model.layer1.bv))

# response to each visible unit being active at 10 and all other
# units being 0. Consistent with normalized spectrogram for pure tone.
# heatmap(sig.(20model.layer1.w .+ model.layer1.b'))

# x = run_spect(model,tone(1kHz,1s,8000));
# x = reshapefor(x,model.layer1);
# standardize!(x,2);

#=
scale = Float64[]
freqs = 10.^(linspace(2,3.5,20))
for f in freqs
  scale = [scale; tone(f*Hz,1s,8000)]
end
for f in reverse(freqs)
  scale = [scale; tone(f*Hz,1s,8000)]
end

x = run_spect(model,scale)
x = reshapefor(x,model.layer1);
standardize!(x,2);

R = sig.(x*model.layer1.w .+ model.layer1.b')'
heatmap(R)

function learn(x)
  n = length(x)
  μ = sum(x) / (n+1)
  σ = sqrt((1 + sum((x .- μ).^2)) / (n+1))

  Normal(μ,σ)
end

p(x) = pdf(learn(x),x)
H(x) = -sum(p(x).*log.(p(x)))

byentropy = mapslices(H,R,2)
indices = sort(1:length(byentropy),by = i -> byentropy[i])

heatmap(R[indices,:])
heatmap(R[indices[end-60:end],:])
=#

# basic steps - present increasing frequency and bandwidth stimuli
# find some way to characterize responses based on this
# highest response doesn't make sense, since some units are default-on
# instead of default-off units

# perhaps we could calculate the mutual information for each
# unit with a on-off pattern for a given frequency
# then we would find an assignment of units to a frequency bin
# such that we maximized the mutual information, and such that
# each frequency has at least one unit.

# we can then repeat this process only for sounds of increasing
# bandwidth, account for a given unit's "best frequency".



# basic steps for mutual information:
# calculate unit responses to silence and the tone
# calculate response distribution (assuming normal)
# of units given tone presences and tone absence

# TODO: in progress - calculate frequency bin MI
# once that's working, we need to find
# the "best frequency" for each unit.

# TODO: I definitely think it is worth trying to use the sigmoid
# transform in the model...

function learn(x)
  n = length(x)
  μ = sum(x) / (n+1)
  sum((x .- μ).^2)
  σ = sqrt((1 + sum((x .- μ).^2)) / (n+1))

  Normal(μ,σ)
end

freqs = (10.^linspace(2,3.5,40))*Hz
freq_resp = zeros(length(freqs),350);

for (i,f) in enumerate(freqs)
  cur_tone = [silence(1s); tone(f, 1s)]
  l1 = run(model,1,cur_tone,upto=1)

  # mutual information calculation assumes first half is silence and second the
  # tone in question

  # NOTE: it's quite possible I can just check the bias and
  # use that to invert values.

  freq_resp[i,:] = mapslices(l1,1) do response
    half = div(size(response,1),2)

    p_silence = learn(response[1:half])
    p_tone = learn(response[half+1:end])
    p_all = learn(response)

    joint_tone = 0.5pdf(p_tone,response)
    joint_silence = 0.5pdf(p_silence,response)

    sum(joint_tone .* log.(joint_tone ./ 0.5pdf.(p_all,response))) +
      sum(joint_silence .* log.(joint_silence ./ 0.5pdf.(p_all,response)))
  end
end

norm(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x))

freq_weight = mapslices(freq_resp,1) do freqs
  freqs = norm(freqs)
  x = sum((1:length(freqs)) .* freqs) ./ sum(freqs)
  isnan(x) ? 0 : x
end

freq_width = mapslices(norm(freq_resp),1) do freqs
  h = max(0.25,maximum(freqs))
  sum(freqs .> 0.2h)
end

#freq_sort = sort(1:size(freq_resp,2),by=i -> 20freq_width[i])
freq_sort = sort(1:size(freq_resp,2),by=i -> 50freq_width[i] + freq_weight[i])

heatmap(freq_resp[:,freq_sort])
plot!(freq_width[freq_sort])


# okay, now I have a good way to organize layer 1, now show the responses,
# accounting for default off and default on (by looking at biases) and correct
# so all are default off for display purposes (by inverting the final result),
# then show the response of layer 1 to the aba stimulus

seq = aba(60ms,10,1kHz,4)
spect = run(model,1,seq,upto=0)
l1_resp = run(model,1,seq,upto=1)
l1_resp .= ifelse.((model.layer1.b .<= 0)',l1_resp,1 .- l1_resp)

heatmap(l1_resp[:,reverse(freq_sort)]',subplot=1,layout=(2,1))
plot!(model.spect,spect,subplot=2)

# what to do for layer 2?
#
# I want to look at the response of elements to
# single tones, rising and falling tones
# and figure out if there is a way to organize them

scale = Float64[]
freqs = 10.^(linspace(2,3.5,20))
for f in freqs
  scale = [scale; tone(f*Hz,1s)]
end
for f in reverse(freqs)
  scale = [scale; tone(f*Hz,1s)]
end

l2_resps = run(model,1:3,scale,upto=2)

heatmap(l2_resps[1]')
heatmap(l2_resps[2]')
heatmap(l2_resps[3]')
