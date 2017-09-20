using Distributions
using Plots; plotlyjs()

import Unitful: ms, Hz

include("model.jl")
include("stim.jl")
include("audio_spect.jl")

data_dir = "/Users/davidlittle/Data"
fs = 8000

deltas = [1,3,6,9,15]
freqs = [500Hz,750Hz,1000Hz,1250Hz,1500Hz,1750Hz,2000Hz]

ab_repeats = 17
tone_len = 60ms

taus = [1,3];

if !isdefined(:model)
  model = Model("/Users/davidlittle/Data",l1_c_a = 0.0f0)
end

a_dist = zeros(length(deltas),maximum(taus),length(freqs),2)
b_dist = zeros(length(deltas),maximum(taus),length(freqs),2)

@time for (i,freq) in enumerate(freqs)
  for (j,delta) in enumerate(deltas)
    @show (freq,delta)
    for k in 1:2
      if k == 1
        base, a_ref, b_ref = alter_ab(tone_len,ab_repeats,fs,freq,delta)
      else
        base, a_ref, b_ref = sync_ab(tone_len,ab_repeats,fs,freq,delta)
      end

      hebb_base = run(model,taus,base)
      hebb_a_ref = run(model,taus,a_ref)
      hebb_b_ref = run(model,taus,b_ref)

      for tau in taus
        a_dist[j,tau,i,k] =
          vecnorm(hebb_base[tau][end-5:end,:].-hebb_a_ref[tau][end-5:end,:])
        b_dist[j,tau,i,k] =
          vecnorm(hebb_base[tau][end-5:end,:].-hebb_b_ref[tau][end-5:end,:])
      end
    end
  end
end

const norm = Normal()
dprime(a,b) = quantile(norm,a) - quantile(norm,b)

function ab_resp(a_dist,b_dist,taus)
  # asum = sum(mean(a_dist,1)[:,taus,:],2)[:]
  # bsum = sum(mean(b_dist,1)[:,taus,:],2)[:]
  asum = sum(mean(a_dist,3)[:,taus],2)[:]
  bsum = sum(mean(b_dist,3)[:,taus],2)[:]
  t = 0.8
  sigma = 0.08
  N = 100

  a_resp = sum(asum' .< t+sigma*randn(11,length(deltas),N),3)
  b_resp = sum(bsum' .< t+sigma*randn(11,length(deltas),N),3)

  dprime.(a_resp./(N+1),b_resp./(N+1))
end

alter = mean(ab_resp(a_dist[:,:,:,1],b_dist[:,:,:,1],taus),1)[:]
sync = mean(ab_resp(a_dist[:,:,:,2],b_dist[:,:,:,2],taus),1)[:]

plot(deltas,[alter sync],label = ["alternating" "synchronous"])
