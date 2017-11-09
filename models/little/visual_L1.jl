################################################################################
# Setup

import Unitful: ms, Hz, kHz, s
using Plots; plotlyjs()
using Distributions
using InformationMeasures
using RCall
using DataFrames

include("model.jl")
include("stim.jl")

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

N = 40
freqs = (10.^linspace(2,3.5,N))*Hz
sil_resp = zeros(length(freqs),350);
freq_resp = zeros(length(freqs),350);

for (i,f) in enumerate(freqs)
  cur_tone = [silence(1s); tone(f, 1s)]
  l1 = run(model,1,cur_tone,upto=1)

  half = div(size(l1,1),2)
  sil_resp[i,:] = mean(l1[1:half,:],1)
  freq_resp[i,:] = mean(l1[half+1:end,:],1)
end

norm(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x))

diffs = abs.(freq_resp .- sil_resp)
freq_stats = mapslices(diffs,1) do freqs
  n = length(freqs)
  freqs = norm(freqs)
  x = sum((1:length(freqs)) .* freqs) ./ sum(freqs)

  center = isnan(x) ? 0 : x
  center_val = clamp(freqs[clamp(round(Int,center),1,n)],0.2,1)

  bandwidth = sum(freqs .> 0.5center_val)
  if bandwidth == 0
    bandwidth = n
  end
  bandwidth = div(bandwidth,5)

  [center,bandwidth]
end

#freq_sort = sort(1:size(diffs,2),by=i -> 20freq_width[i])
freq_sort = sort(1:size(diffs,2),by=i -> 50freq_stats[2,i] + freq_stats[1,i])

df = @mat_asdf(diffs[:,freq_sort],freq_bin,unit,response)

df2 = DataFrame(bandwidth = freq_stats[2,freq_sort],
                freq = freq_stats[1,freq_sort],
                index=1:length(freq_sort))

dir = "../../plots/meeting_2017_11_09/"

R"""
library(ggplot2)
ggplot($df,aes(x = unit,y = freq_bin,fill = response)) + geom_raster() +
  theme_classic() + xlab('Hidden Unit Index') + ylab('Frequency Bin') +
  scale_fill_gradient2(low='black',mid='red',high='yellow',
    name='|r_silence - r_tone|',midpoint=0.5) +
  geom_line(data=$df2,aes(x=index,y=bandwidth,fill=0),color='black',size=3) +
  geom_line(data=$df2,aes(x=index,y=bandwidth,fill=0),color='white') +
  ggtitle("L1 Unit Organization by Bandwidth and Frequency")
ggsave(paste($dir,"l1_tone_map.pdf",sep=""),width=8,height=6)
"""


for f in [250Hz, 500Hz, 1kHz, 2kHz]
  seq = aba(60ms,10,f,6)
  spect = run(model,1,seq,upto=0)
  l1_resp = run(model,1,seq,upto=1)
  l1_diff = abs.(l1_resp .- sil_resp[1,:]')

  df_l1 = @mat_asdf(l1_diff[:,freq_sort]',bin,time,response)
  df_spect = @mat_asdf(spect,time,bin,response)

  freq_str = string(round(Int,ustrip(uconvert(Hz,f))))*"Hz"
  dir = "../../plots/meeting_2017_11_09/"

R"""
  library(ggplot2)
  library(cowplot)

  g_l1 = ggplot($df_l1,aes(x=time,y=bin,fill=response)) + geom_raster() +
    xlab('Time Bin') + ylab('Unit') +
    scale_fill_gradient2(low='black',mid='red',high='yellow',
      name='|r_x - r_silence|',midpoint=0.5,limits=c(0,1))

  g_spect = ggplot($df_spect,aes(x=time,y=bin,fill=response)) + geom_raster() +
    xlab('Time Bin') + ylab('Frequency Bin') +
    scale_fill_gradient2(low='black',mid='red',high='yellow',
      name='Response',midpoint=10)
  body = plot_grid(g_l1,g_spect,ncol=1,align='v',rel_heights=c(3,1))
  title = ggdraw() + draw_label(paste("L1 ABA (",$freq_str,",df = 6)",sep=""),
                                fontface='bold')

  final = plot_grid(title,body,ncol=1,rel_heights=c(0.1,1))
  save_plot(paste($dir,"l1_aba_",$freq_str,"_delta6.pdf",sep=""),final,
            base_height=6)
"""

end


for Δf in [1,3,6,9]
  seq = aba(60ms,10,500Hz,Δf)
  spect = run(model,1,seq,upto=0)
  l1_resp = run(model,1,seq,upto=1)
  l1_diff = abs.(l1_resp .- sil_resp[1,:]')

  df_l1 = @mat_asdf(l1_diff[:,freq_sort]',bin,time,response)
  df_spect = @mat_asdf(spect,time,bin,response)

  delta_str = string(Δf)
  dir = "../../plots/meeting_2017_11_09/"
R"""
  library(ggplot2)
  library(cowplot)

  g_l1 = ggplot($df_l1,aes(x=time,y=bin,fill=response)) + geom_raster() +
    xlab('Time Bin') + ylab('Unit') +
    scale_fill_gradient2(low='black',mid='red',high='yellow',
      name='|r_x - r_silence|',midpoint=0.5,limits=c(0,1))

  g_spect = ggplot($df_spect,aes(x=time,y=bin,fill=response)) + geom_raster() +
    xlab('Time Bin') + ylab('Frequency Bin') +
    scale_fill_gradient2(low='black',mid='red',high='yellow',
      name='Response',midpoint=10)
  body = plot_grid(g_l1,g_spect,ncol=1,align='v',rel_heights=c(3,1))

  title = ggdraw() +
    draw_label(paste("L1 ABA (1000Hz, df = ",$delta_str,")",sep=""),
               fontface='bold')

  final = plot_grid(title,body,ncol=1,rel_heights=c(0.1,1))
  save_plot(paste($dir,"l1_aba_1000Hz_delta",$delta_str,".pdf",sep=""),final,
            base_height=6)
"""
end
