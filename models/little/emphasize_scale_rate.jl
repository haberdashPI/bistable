push!(LOAD_PATH,"packages")
using ProgressMeter
using AuditoryModel
using AuditoryCoherence
using RCall
include("util/stim.jl")

R"library(cowplot)"

quartz() = R"quartz()"

sparams = Dict(:len=>25)
cparams = Dict(:scales => cycoct.*round.(2.0.^linspace(-1,2,9),1),
               :rates => default_rates)

dir = "../../plots/run_2018_02_28"
isdir(dir) || mkdir(dir)

x = @> ab(120ms,120ms,1,6,500Hz,6) normpower amplify(-20)
sp = audiospect(x;sparams...);
cr = cortical(sp;cparams...);

cohere = CoherenceModel(AuditoryModel.Params(cr),4,window=100ms,method=:nmf,
                        delta=50ms,maxiter=200,tol=1e-3)

################################################################################
# how does weighting the rates go?

# it appears that, almost irrespective of the center rate,
# the rate sigma determines stimulus fusion or streaming

function rate_weight(cr,center,sigma)
  r_weights = exp.(.-(log.(abs.(ustrip.(rates(cr)))) .-
                      log.(abs.(ustrip.(center)))).^2 ./
                   sigma^2*log(2))
  cr .* reshape(r_weights,1,:,1,1)
end

crs = Dict()
@showprogress for r in unique(abs.(rates(cr)))
  crs[r] = rate_weight(cr,r,2)
end

# Cs = Dict()
p = []
for r in unique(abs.(rates(cr)))
  if !haskey(Cs,r); Cs[r] = cohere(crs[r]) end
  push!(p,rplot(cohere,Cs[r],components=1:2,scales=[0.5cycoct,1cycoct,3cycoct]))
end

R"""
ps = $p
rates = $(round.(ustrip.(unique(abs.(rates(cr)))),1))

for(i in 1:length(ps)){
  ps[[i]] = ps[[i]] + ggtitle(paste(rates[i]," Hz Emphasis"))
}

p = do.call(plot_grid,c(ps,list(nrow=3,ncol=3)))
save_plot($(joinpath(dir,"2_fake_rate_emphasis.pdf")),p,
          base_aspect_ratio=1.3,nrow=2,ncol=2)
"""

x = @> ab(60ms,60ms,1,12,500Hz,6) normpower amplify(-20)
sp = audiospect(x;sparams...);
cr = cortical(sp;cparams...);

crs = Dict()
@showprogress for r in unique(abs.(rates(cr)))
  crs[r] = rate_weight(cr,r,2)
end

# Cs2 = Dict()
p = []
for r in unique(abs.(rates(cr)))
  if !haskey(Cs2,r); Cs2[r] = cohere(crs[r]) end
  push!(p,rplot(cohere,Cs2[r],components=1:2,scales=[0.5cycoct,1cycoct,3cycoct]))
end


R"""
ps = $p
rates = $(round.(ustrip.(unique(abs.(rates(cr)))),1))

for(i in 1:length(ps)){
  ps[[i]] = ps[[i]] + ggtitle(paste(rates[i]," Hz Emphasis"))
}

p = do.call(plot_grid,c(ps,list(nrow=3,ncol=3)))
save_plot($(joinpath(dir,"3_fake_rate_emphasis_fast_stimulus.pdf")),p,
          base_aspect_ratio=1.3,nrow=2,ncol=2)
"""

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
