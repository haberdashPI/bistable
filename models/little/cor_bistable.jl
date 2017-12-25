include("units.jl")
include("tempc.jl")
include("stim.jl")
include("adaptmi.jl")
setup_sound(sample_rate=8kHz)

dir = "../../plots/run_2017_12_25"
mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=10)
cort = CorticalModel(spect,scales=2.0.^(-1:4),rates=[-2.0.^(1:5); 2.0.^(1:5)])

x = @>(ab(120ms,120ms,1,10,500Hz,6),attenuate(10))

## adaptation MI in spectrogram

sx = spect(x);

sxa = similar(sx)
k = 5
w = [1-exp(-(i - j)^2/k) for i in 1:size(sx,2), j in 1:size(sx,2)]
w ./= sum(w,1)
W_m(x) = w*x

params = AdaptMI(c_m=5,τ_m=30ms,W_m=W_m,c_a=1,τ_a=1s,
                 shape_y = x -> max(0,x))
sxa,a,m = adaptmi(sxa,params) do y_t,t,dt_y
  sx[t,:]
end

p = rplot(spect,sxa)
R"""
library(cowplot)
save_plot($(dir*"/spect_adapt1.pdf"),$p,base_aspect_ratio=1.5)
"""


p = rplot(spect,a)
R"""
library(cowplot)
save_plot($(dir*"/spect_adapt1_a1.pdf"),$p,base_aspect_ratio=1.5)
"""


p = rplot(spect,m)
R"""
library(cowplot)
save_plot($(dir*"/spect_adapt1_m1.pdf"),$p,base_aspect_ratio=1.5)
"""

params = AdaptMI(c_m=5,τ_m=200ms,W_m=W_m,c_a=2,τ_a=1s,
                 shape_y = x -> max(0,x))
sxa,a,m = adaptmi(sxa,params) do y_t,t,dt_y
  sx[t,:]
end

p = rplot(spect,sxa)
R"""
library(cowplot)
save_plot($(dir*"/spect_adapt2.pdf"),$p,base_aspect_ratio=1.5)
"""


p = rplot(spect,a)
R"""
library(cowplot)
save_plot($(dir*"/spect_adapt2_a1.pdf"),$p,base_aspect_ratio=1.5)
"""


p = rplot(spect,m)
R"""
library(cowplot)
save_plot($(dir*"/spect_adapt2_m1.pdf"),$p,base_aspect_ratio=1.5)
"""

########################################
# how do these changes in the spectral layer affect the separation output

# TODO: make better variable names
tempc = TCAnalysis(cort,1,1s)
cr = [cort(sxa),cort(sx)]
Csxa = [tempc(cri) for cri in cr]
signal = [fusion_signal(tempc,Csxa_i,cri) for (cri,Csxa_i) in zip(cr,Csxa)]

# TODO: define this function in seperate script, not in cort_test.jl

p = plot_resps(signal,["yes", "no"],"adaptation/MI")
R"""
library(cowplot)
save_plot($(dir*"/adapt_separation.pdf"),
  $p + ggtitle("Effects of Adaptation and MI (for df = 6st, dt = 240ms)"),
  base_aspect_ratio=1.5)
"""

# what do the components look like
p = Array{Any}(2)
p[1] = rplot(tempc,Csxa[1][round(Int,4s / Δt(spect))])
p[2] = rplot(tempc,Csxa[2][round(Int,4s / Δt(spect))])


R"""
library(cowplot)

p = plot_grid($(p[1]) + ggtitle("Component with Adaptation & MI (4 s)"),
          $(p[2]) + ggtitle("Component without Adaptation & MI (4 s)"),
          align = "h", ncol = 2)
save_plot($(dir*"/adapt_components.pdf"),p,base_aspect_ratio=1.3,ncol=2,nrow=1)
"""
