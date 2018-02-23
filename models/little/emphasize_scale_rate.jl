push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

quartz() = R"quartz()"

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",
                            len=25,min_freq = 400.0Hz,max_freq = 2kHz)
cort = CorticalModel(spect,scales=2.0.^linspace(-0.5,2,6))
cohere = CoherenceModel(cort,4,window=100ms,method=:nmf,delta=50ms,
                       maxiter=200,tol=1e-3)

dir = "../../plots/run_2018_02_20"
isdir(dir) || mkdir(dir)

x = @> ab(120ms,120ms,1,6,500Hz,6) normpower amplify(-20)
sp = spect(x);
cr = cort(sp);

################################################################################
# how does the scales go?
function scale_weight(cr,center,sigma)
  s_weights = exp.(.-(log.(scales(cort)) .- log.(center)).^2 ./ sigma^2*log(2))
  cr .* reshape(s_weights,1,1,:,1)
end

cr = cort(spect(x));
cs05 = scale_weight(cr,0.71,0.1);
cs1 = scale_weight(cr,4,0.1);

C_s05 = cohere(cs05);
C_s1 = cohere(cs1);

p1 = rplot(cohere,C_s05)
p2 = rplot(cohere,C_s1)

sp1 = mean_spect2(cohere,C_s05,cr,component=1)
sp2 = mean_spect2(cohere,C_s1,cr,component=1)

p3 = rplot(spect,sp1)
p4 = rplot(spect,sp2)

R"""
p = plot_grid($p1 + ggtitle("Components for Artificially Emphasized 0.71 cyc/oct"),
              $p2 + ggtitle("Components for Artificially Emphasized 4 cyc/oct"),
              $p3 + ggtitle("Windowed Masking (0.71 cyc/oct)"),
              $p4 + ggtitle("Windowed Masking (4 cyc/oct)"),
              nrow=4,ncol=1)

save_plot($(joinpath(dir,"2_fake_scale_emphasis.pdf")),p,
  base_aspect_ratio=1.3,nrow=4,ncol=2)
"""

################################################################################
# TODO: update below, if I use it
################################################################################
# how does weighting the rates go?

# it appears that, almost irrespective of the center rate,
# the rate sigma determines stimulus fusion or streaming

function rate_weight(cr,center,sigma)
  r_weights = exp.(.-(log.(abs.(rates(cort))) .- log.(abs.(center))).^2 ./
                   sigma^2*log(2))
  cr .* reshape(r_weights,1,:,1,1)
end

cr = cort(spect(x));
ct32_s1 = rate_weight(cr,32,1);
ct1_s1 = rate_weight(cr,1,1);

ct32_s5 = rate_weight(cr,32,5);
ct1_s5 = rate_weight(cr,1,5);

Ct32_s1 = cohere(ct32_s1);
Ct1_s1 = cohere(ct1_s1);

Ct32_s5 = cohere(ct32_s5);
Ct1_s5 = cohere(ct1_s5);

p1 = rplot(cohere,Ct32_s1)
p2 = rplot(cohere,Ct1_s1)
p3 = rplot(cohere,Ct32_s5)
p4 = rplot(cohere,Ct1_s5)

R"""
p = plot_grid($p1 + ggtitle("32 Hz Emphasis, Sigma = 1"),
              $p2 + ggtitle("1 Hz Emphasis, Sigma = 1"),
              $p3 + ggtitle("32 Hz Emphasis, Sigma = 5"),
              $p4 + ggtitle("1 Hz Emphasis, Sigma = 5"),
              nrow=2,ncol=2)
save_plot($(joinpath(dir,"3_fake_rate_emphasis.pdf")),p,
          base_aspect_ratio=1.3,nrow=2,ncol=2)
"""
