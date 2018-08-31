library(cowplot)
library(moments)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)
source("util.R")

dir = file.path("..","..","..","plots",paste("scale_percept_lengths_",
                                             Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)
df = read_feather(file.path("..","..","..","data","count_lengths",
                            "scale_percept_lengths_2018-07-25.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "params_2018-07-24.feather"))
params$pindex = 1:nrow(params)

cleaned_df = df %>% group_by(pindex,created) %>%
  mutate(is_bound = set_bound(stimulus))

summary = cleaned_df %>% group_by(pindex) %>%
  summarize(num_sims      = length(unique(created)),
            N             = sum(!is_bound),
            W             = W(log10(length[!is_bound])),
            kurt          = kurtosis(log10(length[!is_bound])),
            skewness      = skewness(log10(length[!is_bound])),
            mean_ratio    = clean_ratio(stimulus,length,is_bound))

norm_ratio = function(x){
  span = quantile(x[!is.infinite(x)],c(0.05,0.95),na.rm=T)
  (x / pmax(1e-8,span[2] - span[1],na.rm=T))
}

sim_len = df %>%
  filter(pindex == 1,created == first(created)) %>%
  select(length) %>% sum

summary = summary %>%
  mutate(N_per_sec_n   = norm_ratio(N / num_sims/sim_len),
         mean_ratio_n  = norm_ratio(mean_ratio),
         match         = pmax(W,mean_ratio_n)) %>%
  left_join(params) %>%
  gather(measure,value,N:match)

########################################
# N looks pretty similar across noise levels and
# delta_f
p = ggplot(filter(summary,measure == "N",condition == "scales"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(delta_f~θ,labeller=
             label_bquote(rows = Delta[f] == .(delta_f),
                          cols = paste("Threshold ", (theta == .(θ))))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")
save_plot(file.path(dir,"bistable_scale_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=6,base_width=2,base_height=2)

########################################
# selectivity

p1 = ggplot(filter(summary,measure == "mean_ratio",condition == "scales",
                   c_σ == 0.2),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(delta_f~θ,labeller=
             label_bquote(rows = Delta[f] == .(delta_f),
                          cols = paste("Threshold ", (theta == .(θ))))) +
  scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                       breaks=c(-1,0,1),
                       labels=c("-1 (10x shorter)",
                                " 0 (Equal Lengths)",
                                " 1 (10x longer)")) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle("Ratio of Percept Lengths")
p1

fuse_split = summary %>%
  filter(measure == "mean_ratio", condition == "scales", c_σ == 0.2) %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(delta_f,value) %>%
  mutate(quality = -(abs(`3`) + abs(1-`0.5`) + abs(-1 - `12`)))

p2 = ggplot(fuse_split,
            aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=quality)) +
  geom_raster() +
  facet_grid(~θ,labeller=
             label_bquote(cols = paste("Threshold ", (theta == .(θ))))) +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-3,0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",Delta[f],"= 3: -(",
                           abs(Delta[3]) - abs(1-Delta[0.5]) -
                             abs(-1 - Delta[12]),") ")))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align="v")
save_plot(file.path(dir,"bistable_scales_selectivity.pdf"),p,
          base_height=2,base_width=2,nrow=4,ncol=6)


