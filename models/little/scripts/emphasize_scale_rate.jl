# include("units.jl")
include("tempc.jl")
include("stim.jl")

quartz() = R"quartz()"

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",
                            len=25,min_freq = 400.0Hz,max_freq = 2kHz)
cort = CorticalModel(spect,scales = 2.0.^linspace(-2,1,10))
tempc = TCAnalysis(cort,4,window=600ms,method=(:real_pca,8),frame_len=500ms)

dir = "../../plots/run_2018_02_01"
isdir(dir) || mkdir(dir)

x = @> ab(120ms,120ms,1,10,500Hz,6) normpower amplify(-20)
sp = spect(x);
cr = cort(sp);

################################################################################
# how does the scales go?
function scale_weight(cr,center,sigma)
  s_weights = exp.(.-(log.(scales(cort)) .- log.(center)).^2 ./ sigma^2*log(2))
  cr .* reshape(s_weights,1,1,:,1)
end

cr = cort(spect(x));
cs05 = scale_weight(cr,0.35,0.1);
cs1 = scale_weight(cr,1.59,0.1);

# C = tempc(cr);
C_s05 = tempc(cs05);
C_s1 = tempc(cs1);
sp1 = mean_spect(tempc,C_s05,cr,component=1)
sp2 = mean_spect(tempc,C_s1,cr,component=1)


p1 = rplot(tempc,C_s05[3s],n=2,showvar=false)
p2 = rplot(tempc,C_s1[3s],n=2,showvar=false)
p3 = rplot(spect,sp1)
p4 = rplot(spect,sp2)

R"""
p = plot_grid($p1 + ggtitle("Artificially Emphasized 0.25 cyc/oct (at 3s)"),
              $p2 + ggtitle("Artificially Emphasized 1.25 cyc/oct (at 3s)"),
              $p3 + ggtitle("Windowed Masking (0.25 cyc/oct)"),
              $p4 + ggtitle("Windowed Masking (1.25 cyc/oct)"),
              nrow=4,ncol=1)

save_plot($(joinpath(dir,"1_fake_scale_emphasis.pdf")),p,
  base_aspect_ratio=1.3,nrow=4,ncol=2)
"""

################################################################################
# TODO: update below, if I use it
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
save_plot($(joinpath(dir,"rates_6st.png")),p,base_aspect_ratio=1.5)
"""

p = rplot(tempc,C_r16[3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st_rate16.png")),$p,base_aspect_ratio=1.2,ncol=2)
"""

p = rplot(tempc,C_r2[3s])
R"""
library(cowplot)
save_plot($(joinpath(dir,"component_6st_rate2.png")),$p,base_aspect_ratio=1.2,ncol=2)
"""

p = Array{Any}(2)
masked_r16 = mask(tempc,C_r16[3s],cr,0);
sp = inv(cort,masked_r16,usematlab=true)
p[1] = rplot(spect,sp)

masked_r2 = mask(tempc,C_r2[3s],cr,π);
sp = inv(cort,masked_r2,usematlab=true)
p[2] = rplot(spect,sp)
R"""
library(cowplot)
p = plot_grid($(p[1]) + ggtitle("Favored Rate: 16Hz"),
          $(p[2]) + ggtitle("Favored Rate: 2Hz"),align="h",ncol=2)
save_plot($(joinpath(dir,"masked_rates.png")),p,
  base_aspect_ratio=1.4,ncol=2)
"""
