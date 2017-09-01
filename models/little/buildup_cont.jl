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
duration = 10s
durations = linspace(500.0ms,duration,20)

taus = 1:3

aba_seq_offset = 4tone_len

audio_spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5")
if !isdefined(:stims) ||
  !isdefined(:old_deltas) ||
  !isdefined(:old_freqs) ||
  !isdefined(:old_durations) ||
  any(old_deltas .!= deltas) ||
  any(old_freqs .!= old_freqs) ||
  any(durations .!= old_durations)

  global old_durations = durations
  global old_deltas = deltas
  global old_freqs = freqs

  info("Initializing Stimuli")
  stims = Array{Any}(length(freqs),length(deltas))

  @time for (i,freq) in enumerate(freqs)
    for (j,delta) in enumerate(deltas)
      @show (freq,delta)
      repeats = ceil(Int,duration / 240ms)+1
      stim = aba(tone_len,repeats,fs,freq,delta)

      stims[i,j] = run(audio_spect,stim)
    end
  end
end

hebb_dist = zeros(length(freqs),length(deltas),length(durations),
                  maximum(taus),2)

function time_to_frame(time,tau)
  frame_len = s/fs * frame_length(audio_spect) * frame_length(model) * tau
  floor(Int,time / frame_len)
end

# TODO: try looking at the tau's sperately, see if
# we should do something different wiht the taus

info("Simulating Streaming...")
@time for (i,freq) in enumerate(freqs)
  for (j,delta) in enumerate(deltas)
    @show (freq,delta)

    hebb = run(model,taus,stims[i,j])
    for (k,dur) in enumerate(durations)
      # compare most recent a at the start of a triplet to the most recent b
      for h in taus
        N_triplets = floor((dur-aba_seq_offset)/(4tone_len))
        a_onset = aba_seq_offset + 4tone_len*N_triplets

        ap_onset = a_onset + 2tone_len
        ap_onset = ap_onset > dur ? ap_onset - 4tone_len : ap_onset

        b_onset = a_onset + tone_len
        b_onset = b_onset > dur ? b_onset - 4tone_len : b_onset

        afrom = time_to_frame(a_onset,h)
        apfrom = time_to_frame(ap_onset,h)
        bfrom = time_to_frame(b_onset,h)

        len = max(1,time_to_frame(tone_len,h))

        a_range = afrom:(afrom+len-1)
        ap_range = apfrom:(apfrom+len-1)
        b_range = bfrom:(bfrom+len-1)

        hebb_dist[i,j,k,h,:] =
          [vecnorm(hebb[h][a_range,:]-hebb[h][b_range,:]),
           vecnorm(hebb[h][a_range,:]-hebb[h][ap_range,:])]
      end
    end
  end
end

normed = hebb_dist[:,:,:,1] ./ hebb_dist[:,:,:,2]
mdist = squeeze(mean(normed,1),1)

responses = respond(mdist,1.6,0.4^2)

plot(ustrip(durations),responses',label=deltas')
