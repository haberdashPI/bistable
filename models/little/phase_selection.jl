push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
using DataFrames
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"
dir = "../../plots/run_2018_02_01"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=10,
                            min_freq = 250Hz,max_freq=1500Hz)
cort = CorticalModel(spect,rates=sort([-2.^(-1:0.5:5); 2.^(-1:0.5:5)]))
tempc = TCAnalysis(cort,4,window=600ms,method=(:real_pca,8),frame_len=500ms)

x_a = ab(120ms,120ms,1,10,500Hz,6,:without_a) |> normpower |> amplify(-20dB)
x_b = ab(120ms,120ms,1,10,500Hz,6,:without_b) |> normpower |> amplify(-20dB)
x = ab(120ms,120ms,1,10,500Hz,6) |> normpower |> amplify(-20dB)

sp = spect(x);
cr = cort(sp);
C = tempc(cr);

# how do the real-valued components work
p1 = rplot(tempc,C[3s],λ_digits=3,showvar=false,n=3)

sp_C = mean_spect(tempc,C,cr,component=1)
p2 = rplot(spect,sp_C)

sp_C = mean_spect(tempc,C,cr,component=2)
p3 = rplot(spect,sp_C)

sp_C = mean_spect(tempc,C,cr,component=3)
p4 = rplot(spect,sp_C)


R"""
p = plot_grid($p1 + ggtitle(expression(paste(Delta*f==6," masks at 3 seconds"))),
          $p2 + ggtitle("Windowed Masking (1st component)"),
          $p3 + ggtitle("Windowed Masking (2nd component)"),
          $p4 + ggtitle("Windowed Masking (3rd component)"),
          nrow=4)
save_plot($(joinpath(dir,"1_real_masks_df6.pdf")),p,
  base_aspect_ratio=2,nrow=4)
"""

# how do they respond to varying Δf
sp12 = ab(120ms,120ms,1,10,500Hz,12) |> normpower |> amplify(-20dB) |> spect
C12 = tempc(sp12)
p1 = rplot(tempc,C12[3s],n=3,showvar=false)

sp_C = mean_spect(tempc,C12,cort(sp12),component=1)
p2 = rplot(spect,sp_C)

sp_C = mean_spect(tempc,C12,cort(sp12),component=2)
p3 = rplot(spect,sp_C)

sp_C = mean_spect(tempc,C12,cort(sp12),component=3)
p4 = rplot(spect,sp_C)


R"""
p = plot_grid($p1 + ggtitle(expression(paste(Delta*f==12," masks at 3 seconds"))),
          $p2 + ggtitle("Windowed Masking (1st component)"),
          $p3 + ggtitle("Windowed Masking (2nd component)"),
          $p4 + ggtitle("Windowed Masking (3rd component)"),
          nrow=4)
save_plot($(joinpath(dir,"2_real_masks_df12.pdf")),p,
  base_aspect_ratio=2,nrow=4)
"""

sp1 = ab(120ms,120ms,1,10,500Hz,1) |> normpower |> amplify(-20dB) |> spect
C1 = tempc(sp1)
p1 = rplot(tempc,C1[3s],n=3,showvar=false)

sp_C = mean_spect(tempc,C1,cort(sp1),component=1)
p2 = rplot(spect,sp_C)

sp_C = mean_spect(tempc,C1,cort(sp1),component=2)
p3 = rplot(spect,sp_C)

sp_C = mean_spect(tempc,C1,cort(sp1),component=3)
p4 = rplot(spect,sp_C)


R"""
p = plot_grid($p1 + ggtitle(expression(paste(Delta*f==1," masks at 3 seconds"))),
          $p2 + ggtitle("Windowed Masking (1st component)"),
          $p3 + ggtitle("Windowed Masking (2nd component)"),
          $p4 + ggtitle("Windowed Masking (3rd component)"),
          nrow=4)
save_plot($(joinpath(dir,"3_real_masks_df1.pdf")),p,
  base_aspect_ratio=2,nrow=4)
"""

# how do they respond to varying Δd
sp60 = ab(60ms,60ms,1,20,500Hz,6) |> normpower |> amplify(-20dB) |> spect
C60 = tempc(sp60)
p1 = rplot(tempc,C60[3s],n=3,showvar=false)

p2 = rplot(tempc,C[3s],λ_digits=3,showvar=false,n=3)

sp240 = ab(240ms,240ms,1,5,500Hz,6) |> normpower |> amplify(-20dB) |> spect
C240 = tempc(sp240)
p3 = rplot(tempc,C240[3s],n=3,showvar=false)

R"""
p = plot_grid($p1 + ggtitle(expression(paste(Delta*t==60,"ms masks at 3 seconds"))),
              $p2 + ggtitle(expression(paste(Delta*t==120,"ms masks at 3 seconds"))),
              $p3 + ggtitle(expression(paste(Delta*t==240,"ms masks at 3 seconds"))),
          nrow=3)
save_plot($(joinpath(dir,"4_real_masks_dt_all.pdf")),p,
  base_aspect_ratio=2,nrow=3)
"""

# emphasize particular scales
function scale_weight(cr,center)
  s_weights = exp.(.-(log.(scales(cort)) .- log.(center)).^2 ./ 0.1log(2))
  cr .* reshape(s_weights,1,1,:,1)
end
cs1 = scale_weight(cr,1);
cs025 = scale_weight(cr,0.25);

Cs1 = tempc(cs1)
p1 = rplot(tempc,Cs1[3s],n=2,showvar=false)
sp_Cs1 = mean_spect(tempc,Cs1,cr,component=1)
p2 = rplot(spect,sp_Cs1)

Cs025 = tempc(cs025)
p3 = rplot(tempc,Cs025[3s],n=2,showvar=false)
sp_Cs025 = mean_spect(tempc,Cs025,cr,component=1)
p4 = rplot(spect,sp_Cs025)

R"""
p = plot_grid($p1 + ggtitle("Scale 1 cyc/oct Component"),
              $p2 + ggtitle("Scale 1 cyc/oct windowed masking"),
              $p3 + ggtitle("Scale 0.25 cyc/oct Component"),
              $p4 + ggtitle("Scale 0.25 cyc/oct windowed masking"),
              nrow=2,ncol=2)
# save_plot($(joinpath(dir,"5_real_masks_scale_weighting.pdf")),p,
  # base_aspect_ratio=1.6,nrow=2,ncol=2)
"""

function rate_weight(cr,center,σ=0.1)
    r_weights = exp.(.-(log.(abs.(rates(cort))) .- log.(abs.(center))).^2 ./
                     (σ^2*log(2)))
  cr .* reshape(r_weights,1,:,1,1)
end
cr32 = rate_weight(cr,32);
cr1 = rate_weight(cr,0.5);

Cr32 = tempc(cr32)
p1 = rplot(tempc,Cr32[3s],n=2,showvar=false)
sp_Cr32 = mean_spect(tempc,Cr32,cr,component=2)
p2 = rplot(spect,sp_Cr32)

Cr1 = tempc(cr1)
p3 = rplot(tempc,Cr1[3s],n=2,showvar=false)
sp_Cr1 = mean_spect(tempc,Cr1,cr,component=2)
p4 = rplot(spect,sp_Cr1)

R"""
p = plot_grid($p1 + ggtitle("Rate 32 Hz Components"),
              $p2 + ggtitle("Rate 32 Hz windowed masking (component 2)"),
              $p3 + ggtitle("Rate 0.5 Hz Components"),
              $p4 + ggtitle("Rate 0.5 Hz windowed masking (component 2)"),
              nrow=2,ncol=2)
save_plot($(joinpath(dir,"6_real_masks_rate_weighting.pdf")),p,
  base_aspect_ratio=1.6,nrow=2,ncol=2)
"""
