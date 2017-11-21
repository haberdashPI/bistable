import Unitful: ms, Hz, kHz, s
using Plots; plotlyjs()
using Parameters
using DataFrames
using Feather
using RCall

include("model.jl")
include("stim.jl")

# TODO: figure out bug that means after we change l3 params
# we can't go back to baseline for l3

# model progression,
#
# 1. all use_sig = false (with and without l3 threshold function change)
# 2. use_sig l1 = true, l2,l3 = false
# 3. use_sig l1,l2 = true, l3 = false, l3thresh = 0.4
# 4. use_sig l1,l2 = true, l3 = false, l3thresh = 0.1

# TODO: create a function that run this benchmark
# rather than running this as a script

# if !isdefined(:model)
model = Model("/Users/davidlittle/Data",
              l1params = LayerParams(c_mi = 0,c_a = 0,use_sig=true),
              l2params = LayerParams(c_mi = 0,c_a = 0,use_sig=false),
              l3params = LayerParams(c_mi = 0,c_a = 0,use_sig=false))
# end

tone_len = 60ms
deltas = [1, 3, 6, 9]
# deltas = [7, 8, 8.5, 9, 10, 11]
freqs = [1000Hz]#[500Hz 750Hz 1000Hz 1250Hz 1500Hz 1750Hz 2000Hz]
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
    hebb = run(model,taus,aba(tone_len,repeats,freq,delta))
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


# model = Model("/Users/davidlittle/Data",
#               l1params = LayerParams(c_mi = 0,c_a = 0),
#               l2params = LayerParams(c_mi = 0,c_a = 0))
# l1,l2,l3 = run(model,taus,aba(tone_len,32,1000Hz,9),
#                return_all=true)



Base.squeeze(f,A,dims) = squeeze(f(A,dims),dims)

shebb = squeeze(sum,hebb_dist[:,:,:,:,:],4)
normed = shebb[:,:,:,1] ./ shebb[:,:,:,2]
responses = respond(normed,1.0,0.4^2)
plot(ustrip(durations),responses[1,:,:]',label=deltas')

plot(ustrip(durations),normed[1,:,:]',label=deltas')

r_indices = CartesianRange(size(normed))
df = DataFrame(response = normed[:],
               freq = ustrip([freqs[ii[1]] for ii in r_indices][:]),
               delta = [deltas[ii[2]] for ii in r_indices][:],
               time = ustrip([durations[ii[3]] for ii in r_indices][:]))

R"""
library(ggplot2)

ggplot($df,aes(x=time,y=response,color=factor(delta))) +
  geom_line(size=3,color='black',aes(group=delta)) + geom_line(size=1.5) +
  xlab('duration (s)') + ylab('normed dist') +
  theme_classic(base_size = 23) +
  scale_color_brewer(palette='YlGn',name='Delta F (st)') +
  theme()

# ggsave(paste('../../plots/buildup_D_1Hz_L1a_mi_L2a_mi_L3a_mi_',Sys.Date(),'.pdf',sep=''),
#     width=8,height=6)
"""

# plot(ustrip(durations),normed[1,:,:]',label=deltas')

# plot(ustrip(durations),shebb[1,:,:,1]',label=deltas')
# plot(ustrip(durations),shebb[1,:,:,2]',label=deltas')

# responses = respond(normed,1.8,0.5^2)
# plot(ustrip(durations),responses_adapt[3,:,:]',
#      label=deltas',subplot=1,layout=(2,1))

# plot(ustrip(durations),responses_all[3,:,:]',label=deltas',subplot=1,layout=(3,1))
# plot!(ustrip(durations),responses_adapt[3,:,:]',label=deltas',subplot=2,layout=(3,1))
# plot!(ustrip(durations),responses_none[3,:,:]',label=deltas',subplot=3)

# # plot(ustrip(durations),squeeze(mean,responses,1)',label=deltas')

# r_indices = CartesianRange(size(responses))
# df = DataFrame(response = responses[:],
#                freq = ustrip([freqs[ii[1]] for ii in r_indices][:]),
#                delta = [deltas[ii[2]] for ii in r_indices][:],
#                time = ustrip([durations[ii[3]] for ii in r_indices][:]))

# filename = ("../../data/buildup_cont_adapt_"*
#               Dates.format(now(),"yyyy-mm-dd_HH.MM")*".feather")
# Feather.write(filename,df)
# println("Wrote results to "*filename)