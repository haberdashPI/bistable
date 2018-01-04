include("units.jl")
include("tempc.jl")
include("stim.jl")
setup_sound(sample_rate=8kHz)

quartz() = R"quartz()"

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",
                            len=25,min_freq = 400.0Hz,max_freq = 2kHz)
sts = [0.1,1,6,12]

cort = CorticalModel(spect,scales=2.^(-1:0.5:4),bandonly=false)
dir = "../../plots/run_2017_12_24"
mkdir(dir)

f_ab(st) = @>(ab(120ms,120ms,1,10,500Hz,st),attenuate(10))
ab_f = [f_ab(st) for st in sts]
crs = [cort(spect(ab_f[i]),usematlab=true) for i in eachindex(ab_f)];

tempc = TCAnalysis(cort,1,1s,method=:pca)
Cs = [tempc(crs[i]) for i in eachindex(ab_f)];

p = rplot(tempc,Cs,"delta f (st)" => sts)
R"""
library(cowplot)
p = $p + ylab(expression(lambda[1] / sigma^2)) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1.0))
save_plot($(joinpath(dir,"delta_f.pdf")),p,base_aspect_ratio=1.5)
"""

p = rplot(tempc,Cs[3][3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st.pdf")),$p,base_aspect_ratio=1.2,ncol=2)
"""

masked = mask(tempc,Cs[3][3s],crs[3],0);
sp = inv(cort,masked,usematlab=true)
p = rplot(spect,sp)
R"""
library(cowplot)
save_plot($(joinpath(dir,"masked_phase_0_6st.pdf")),$p,base_aspect_ratio=1.5)
"""

masked = mask(tempc,Cs[3][3s],crs[3],π);
sp = inv(cort,masked,usematlab=true)
p = rplot(spect,sp)
R"""
library(cowplot)
save_plot($(joinpath(dir,"masked_phase_pi_6st.pdf")),$p,base_aspect_ratio=1.5)
"""

################################################################################
# how does weighting the rates go?
function rate_weight(cr,center)
  r_weights = exp.(.-(log.(abs.(rates(cort))) .- log.(abs.(center))).^2 ./ 0.5log(2))
  cr .* reshape(r_weights,1,:,1,1)
end

cr = cort(spect(ab_f[3]),usematlab=true);
cr_r16 = rate_weight(cr,16);
cr_r2 = rate_weight(cr,2);
rplot(cort,cr_r16,rates=[-16,-8,-2,2,8,16],scales=[2,4,16])

C = tempc(cr);
C_r16 = tempc(cr_r16);
C_r2 = tempc(cr_r2);

p = rplot(tempc,[C,C_r16,C_r2],"Favored Rate" => ["none","16 Hz","2 Hz"])
R"""
library(cowplot)
p = $p + ylab(expression(lambda[1] / sigma^2)) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1.0))
save_plot($(joinpath(dir,"rates_6st.pdf")),p,base_aspect_ratio=1.5)
"""

p = rplot(tempc,C_r16[3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st_rate16.pdf")),$p,base_aspect_ratio=1.2,ncol=2)
"""

p = rplot(tempc,C_r2[3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st_rate2.pdf")),$p,base_aspect_ratio=1.2,ncol=2)
"""

p = Array{Any}(2)
masked_r16 = mask(tempc,C_r16[3s],cr,0);
sp = inv(cort,masked_r16,usematlab=true)
p[1] = rplot(spect,sp)

masked_r2 = mask(tempc,C_r2[3s],cr,0);
sp = inv(cort,masked_r2,usematlab=true)
p[2] = rplot(spect,sp)
R"""
library(cowplot)
p = plot_grid($(p[1]) + ggtitle("Favored Rate: 16Hz"),
          $(p[2]) + ggtitle("Favored Rate: 2Hz"),align="h",ncol=2)
save_plot($(joinpath(dir,"masked_rates.pdf")),p,
  base_aspect_ratio=1.4,ncol=2)
"""

################################################################################
# how does the scales go?
function scale_weight(cr,center)
  s_weights = exp.(.-(log.(scales(cort)) .- log.(center)).^2 ./ 0.5log(2))
  cr .* reshape(s_weights,1,1,:,1)
end

cr = cort(spect(ab_f[3]),usematlab=true);
cs16 = scale_weight(cr,16);
cs05 = scale_weight(cr,0.5);

# C = tempc(cr);
C_s16 = tempc(cs16);
C_s05 = tempc(cs05);

p = rplot(tempc,[C,C_s16,C_s05],
          "Favored Scale" => ["none","16 cyc/oct","0.5 cyc/oct"])
R"""
library(cowplot)
p = $p + ylab(expression(lambda[1] / sigma^2)) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1.0))
save_plot($(joinpath(dir,"scales_6st.pdf")),p,base_aspect_ratio=1.5)
"""


p = rplot(tempc,C_s16[3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st_scale16.pdf")),$p,base_aspect_ratio=1.2,ncol=2)
"""

p = rplot(tempc,C_s05[3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st_scale05.pdf")),$p,base_aspect_ratio=1.2,ncol=2)
"""


p = Array{Any}(2)
masked_s16 = mask(tempc,C_s16[3s],cr,0);
sp = inv(cort,masked_s16,usematlab=true)
p[1] = rplot(spect,sp)

masked_s05 = mask(tempc,C_s05[3s],cr,0);
sp = inv(cort,masked_s05,usematlab=true)
p[2] = rplot(spect,sp)

R"""
library(cowplot)
p = plot_grid($(p[1]) + ggtitle("Favored Scale: 16 cyc/oct"),
          $(p[2]) + ggtitle("Favored Scale: 0.5 cyc/oct"),align="h",ncol=2)
save_plot($(joinpath(dir,"masked_scales.pdf")),p,
  base_aspect_ratio=1.4,ncol=2)
"""

# TODO: try examining these separation components when
# applying adaptation and MI to the first layer

# TODO: try examining the effects of the weighting across
# different delta f's, is there a clear point of ambiguity
# in separation that happens


################################################################################
# what if we empahsize 1 cyc/oct?
function scale_weight(cr,center)
  s_weights = exp.(.-(log.(scales(cort)) .- log.(center)).^2 ./ 0.5log(2))
  cr .* reshape(s_weights,1,1,:,1)
end

cr = cort(spect(ab_f[3]),usematlab=true);
cs16 = scale_weight(cr,16);
cs1 = scale_weight(cr,1);

C = tempc(cr);
C_s16 = tempc(cs16);
C_s1 = tempc(cs1);

p = rplot(tempc,[C,C_s16,C_s1],
          "Favored Scale" => ["none","16 cyc/oct","1 cyc/oct"])
R"""
library(cowplot)
p = $p + ylab(expression(lambda[1] / sigma^2)) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1.0))
save_plot($(joinpath(dir,"scales_6st.pdf")),p,base_aspect_ratio=1.5)
"""
