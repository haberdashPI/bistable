push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

R"library(ggplot2)"
R"library(cowplot)"
quartz() = R"quartz()"

dir = "../../plots/run_2018_02_01"
isdir(dir) || mkdir(dir)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect)
tempc = TCAnalysis(cort,4,window=750ms,method=(:real_pca,8),frame_len=50ms)

x = ab(120ms,120ms,1,40,500Hz,6) |> normpower |> amplify(-10dB)

sp = spect(x);
cr = cort(sp);
C = tempc(cr);

Ca = similar(C)

# ########################################
# # can one component win?

# params = AdaptMI(c_m=20,τ_m=200ms,
#                  c_a=0,τ_a=2s,shape_y = x -> max(0,x),
#                  Δt = Δt(tempc));

# Ca,a,m = (adaptmi(Ca,params) do C_t,t,dt_C C[t] end);

# p2 = rplot(tempc,Ca)

########################################
# can we we get cycles?

# params = AdaptMI(α = 1.0,
#                  c_m=10,τ_m=400ms,
#                  c_a=4,τ_a=2s,shape_y = x -> max(0,x),
#                  Δt = Δt(tempc));

# Ca,a,m = (adaptmi(Ca,params) do C_t,t,dt_C C[t] end);

params = AdaptMI(α = 0.5,
                 c_m=10,τ_m=400ms,
                 c_a=4,τ_a=2s,shape_y = x -> max(0,x),
                 Δt = Δt(tempc));

p1 = rplot(tempc,C)

Ca,a,m = (adaptmi(Ca,params) do C_t,t,dt_C C[t] end);

p2 = rplot(tempc,Ca)

spa = mean_spect(tempc,Ca,cr,component=:max)
p3 = rplot(spect,spa)

spc1 = mean_spect(tempc,C,cr,component=1)
p4 = rplot(spect,spc1)

p5 = rplot(spect,spa .- spc)



R"""
p = plot_grid($p1 + ggtitle("Unmodified Component Strengths") +
                    ylab(expression(lambda/sigma^2)),
              $p2 + ggtitle("Adapt/MI Component Strengths") +
                    ylab(expression(lambda/sigma^2)),
              $p5 + ggtitle("Spect diff (adapt/mi - unmodified)"),
              $p3 + ggtitle("Windowed Spectrogram by Adapt/MI components."),
              $p4 + ggtitle("Windowed Spectrogram by unmodified components."),
              nrow=5,ncol=1)

save_plot($(joinpath(dir,"2_bistable_components.pdf")),p,
  base_aspect_ratio=3,nrow=5,ncol=1)
"""
