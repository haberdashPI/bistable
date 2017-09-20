import Unitful: ms, Hz, kHz, s
using Plots; plotlyjs()
using Parameters
using DataFrames
using Feather

include("model.jl")
include("stim.jl")

if !isdefined(:model)
  model = Model("/Users/davidlittle/Data")
end

fs = 8000

tone_len = 60ms
deltas = [1, 3, 6, 9]
freqs = [500Hz 750Hz 1000Hz 1250Hz 1500Hz 1750Hz 2000Hz]
duration = 10s
durations = linspace(500.0ms,duration,20)

# TODO: in the process of looking at different tone lengths
# (which require a different set of tau's, I believe)
taus = 1:2

aba_seq_offset = 4tone_len

hebb_dist = zeros(length(freqs),length(deltas),
                  length(durations),maximum(taus),2)

function time_to_frame(time,tau)
  floor(Int,time / (frame_length(model) * tau))
end

info("Simulating Streaming...")
@time for (i,freq) in enumerate(freqs)
  for (j,delta) in enumerate(deltas)
    @show (freq,delta)

    repeats = ceil(Int,duration / 4tone_len)+1
    hebb = run(model,taus,aba(tone_len,repeats,fs,freq,delta))
    for (k,dur) in enumerate(durations)
      N_triplets = floor((dur-aba_seq_offset)/4tone_len)-1
      if N_triplets < 1 continue end

      a_onset = aba_seq_offset + 4tone_len*N_triplets

      b_onset = a_onset + tone_len
      ap_onset = a_onset + 2tone_len

      # compare most recent a at the start of a triplet to the most recent b
      for tau in taus
        afrom = time_to_frame(a_onset,tau)
        apfrom = time_to_frame(ap_onset,tau)
        bfrom = time_to_frame(b_onset,tau)

        len = max(1,time_to_frame(tone_len,tau))

        a_range = afrom:(afrom+len-1)
        ap_range = apfrom:(apfrom+len-1)
        b_range = bfrom:(bfrom+len-1)

        hebb_dist[i,j,k,tau,:] =
          [vecnorm(hebb[tau][a_range,:]-hebb[tau][b_range,:],1),
           vecnorm(hebb[tau][a_range,:]-hebb[tau][ap_range,:],1)]
      end
    end
  end
end

Base.squeeze(f,A,dims) = squeeze(f(A,dims),dims)

shebb = squeeze(sum,hebb_dist[:,:,:,:,:],4)
normed = shebb[:,:,:,1] ./ shebb[:,:,:,2]
responses = respond(normed,0.9,0.4^2)

# plot(ustrip(durations),normed[3,:,:]',label=deltas')
plot(ustrip(durations),squeeze(mean,responses,1)',label=deltas')

# r_indices = CartesianRange(size(responses))
# df = DataFrame(response = responses[:],
#                freq = ustrip([freqs[ii[1]] for ii in r_indices][:]),
#                delta = [deltas[ii[2]] for ii in r_indices][:],
#                time = ustrip([durations[ii[3]] for ii in r_indices][:]))

# filename = ("../../data/buildup_cont_L1_"*
#               Dates.format(now(),"yyyy-mm-dd_HH.MM")*".feather")
# Feather.write(filename,df)
# println("Wrote results to "*filename)
