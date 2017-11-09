import Unitful: ms, Hz, kHz, s
using Distributions
using InformationMeasures
using RCall
using DataFrames

include("model.jl")
include("stim.jl")

# TODO: look at MI of
# layer2 to each tone, and then to
# rising and falling tones


macro mat_asdf(X,x,y,z)
  quote
    let X = $(esc(X)), ixs = CartesianRange(size(X))
      DataFrame($z = X[:],
                $x = map(x -> x[1],ixs)[:],
                $y = map(x -> x[2],ixs)[:])
    end
  end
end

model = Model("/Users/davidlittle/Data",
              l1params = LayerParams(c_a = 0,c_mi = 0,use_sig=true),
              l2params = LayerParams(c_a = 0,c_mi = 0,use_sig=false),
              l3params = LayerParams(c_a = 0,c_mi = 0,use_sig=false))

################################################################################
# map MI of silence vs. tones

N = 40
freqs = (10.^linspace(2,3.5,N))*Hz
sil_resp = zeros(length(freqs),300);
freq_resp = zeros(length(freqs),300);
MI_resp = zeros(length(freqs),300);

# TODO: this is only for tau = 1
# TODO: we may want to look into
# how this visualization changes when
# we mix the different components
# properly

for (i,f) in enumerate(freqs)
  @show f
  cur_tone = [silence(1s); tone(f, 1s)]
  l2 = run(model,1,cur_tone,upto=2)[1]

  half = div(size(l2,1),2)
  sil_resp[i,:] = mean(l2[1:half,:],1)
  freq_resp[i,:] = mean(l2[half+1:end,:],1)
  MI_resp[i,:] = mapslices(l2,1) do freqs
    label = 1:length(freqs) .> 0.5length(freqs)
    get_mutual_information(freqs,label)
  end
end

norm(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x))
diffs = abs.(freq_resp .- sil_resp)


df = @mat_asdf(diffs,freq_bin,unit,response)

R"""
library(ggplot2)
ggplot($df,aes(x = unit,y = freq_bin,fill = response)) + geom_raster() +
  theme_classic() + xlab('Hidden Unit Index') + ylab('Frequency Bin') +
  scale_fill_gradient2(low='black',mid='red',high='yellow',
    name='|r_silence - r_tone|',midpoint=0.5) +
  ggtitle("L2 tonal responses")
"""

df = @mat_asdf(MI_resp,freq_bin,unit,response)

R"""
library(ggplot2)
ggplot($df,aes(x = unit,y = freq_bin,fill = response)) + geom_raster() +
  theme_classic() + xlab('Hidden Unit Index') + ylab('Frequency Bin') +
  scale_fill_gradient2(low='black',mid='red',high='yellow',
    name='MI (tone vs. silence)',midpoint=0.5) +
  ggtitle("L2 MI")
"""

# TODO: now test responses with rising and falling tones
