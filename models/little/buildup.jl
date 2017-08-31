import Unitful: ms, Hz
using Plots; plotlyjs()

include("model.jl")
include("stim.jl")
include("audio_spect.jl")

audio_spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5")

if !isdefined(:model)
  model = Model("/Users/davidlittle/Data/model.h5")
end

tone_len = 60ms
fs = 8000
deltas = [1, 3, 6, 9];
freqs=[500Hz 750Hz 1000Hz 1250Hz 1500Hz 1750Hz 2000Hz];
durations=[500ms, 1000ms, 1500ms, 2000ms, 2500ms, 3000ms, 3500ms, 4000ms, 4500ms,
           5000ms, 5500ms, 6000ms, 6500ms, 7000ms, 7500ms, 8000ms, 8500ms,
           9000ms, 9500ms, 10000ms];

taus = 1:3

audio_spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5")
if !isdefined(:a_ref)
  a_ref = Array{Any}(length(freqs),length(deltas),length(durations))
  b_ref = Array{Any}(length(freqs),length(deltas),length(durations))

  @time for (i,freq) in enumerate(freqs)
    for (j,delta) in enumerate(deltas)
      for (k,duration) in enumerate(durations)
        @show (freq,delta,duration)
        repeats = ceil(Int,duration / 240ms)
        a, b = buildup_aba(tone_len,repeats,fs,freq,delta)

        a_ref[i,j,k] = run(audio_spect,a)
        b_ref[i,j,k] = run(audio_spect,b)
      end
    end
  end
end

dist(x,y) = sqrt(trace((x-y)*(x-y)'))
hebb_dist = zeros(length(freqs),length(deltas),length(durations))
@time for (i,freq) in enumerate(freqs)
  for (j,delta) in enumerate(deltas)
    for (k,duration) in enumerate(durations)
      @show (freq,delta,duration)

      hebb_a = run(model,taus,a_ref[i,j,k])
      hebb_b = run(model,taus,b_ref[i,j,k])

      for tau in taus
        hebb_dist[i,j,k] +=
          dist(hebb_a[tau][end-5:end,:],hebb_b[tau][end-5:end,:])
      end
    end
  end
end

mdist = mean(hebb_dist,1)
mdist = mdist[:,:,1:1] .- mdist

t = 4
sigma = 1
resps = mdist .> t+sigma*randn(size(mdist)...,10000)
prop_streaming = squeeze(mean(resps,4),(1,4))

plot(ustrip(durations),prop_streaming',label = deltas')
