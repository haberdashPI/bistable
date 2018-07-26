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
                            "scale_percept_lengths_2018-07-21.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "params_2018-07-16.feather"))
params$pindex = 1:nrow(params)

cleaned_df = df %>% group_by(pindex,created) %>%
  mutate(is_bound = set_bound(stimulus))

# TODO: break up by freq and scales

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

################################################################################
# scale plots

########################################
# percepts per simulation:

sim_len = df %>%
  filter(pindex == 1,created == first(created)) %>%
  select(length) %>% sum

# N looks pretty similar across noise levels and
# delta_f
p = ggplot(filter(summary,measure == "N",condition == "scales"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")
save_plot(file.path(dir,"bistable_scale_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=6,base_width=2,base_height=2)

########################################
# good view of the selectivity of bistability

p1 = ggplot(filter(summary,measure == "mean_ratio",condition == "scales"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=
             label_bquote(rows = Delta[f] == .(delta_f),
                          cols = paste("Noise ", (c[sigma] == .(c_σ))))) +
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
  filter(measure == "mean_ratio", condition == "scales") %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(delta_f,value) %>%
  mutate(quality = -(abs(`3`) + abs(1-`0.5`) + abs(-1 - `12`)))

p2 = ggplot(fuse_split,
            aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=quality)) +
  geom_raster() +
  facet_grid(~c_σ,labeller=
             label_bquote(cols = paste("Noise ",(c[sigma]) == .(c_σ)))) +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-3,0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",Delta[f],"= 3: -(",
                           abs(Delta[3]) - abs(1-Delta[0.5]) -
                             abs(-1 - Delta[12]),") ")))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align='v')
save_plot(file.path(dir,"bistable_scales_selectivity.pdf"),p,
          base_height=2,base_width=2,nrow=4,ncol=6)

########################################
# how log-normal are the data?

ggplot(filter(summary,measure == "W",condition == "scales"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Wilks-Shapiro",palette="Greens",direction=1) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = summary %>% filter(measure == "kurt",condition == "scales") %>%
  select(value) %>% abs %>% max(na.rm=T)
ggplot(filter(summary,measure == "kurt"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Kurtosis",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = summary %>% filter(measure == "skewness",condition == "scales") %>%
  select(value) %>% abs %>% max(na.rm=T)
ggplot(filter(summary,measure == "skewness"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Skewness",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

########################################
# let's plot a histogram of some parametesr
sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-18) < 6e-1,abs(c_a-10) < 6e-1,
         delta_f==6,condition == "scales") %>%
  select(pindex) %>% first
selected = df %>% filter(pindex == sindex)
ggplot(selected,aes(x=length)) + geom_histogram(bins=15)
ggplot(selected,aes(x=length)) + geom_density()

# TASKS:

# 1. test/run multi-track bistability
# 2. examine freqs - they are almost somewhat non-seqsical
#    for the extreme stimuli
# 3. run freqs with same params as the scales
# 4. run a narrower range of W_m_σ

################################################################################
# plotting of freq condition

########################################
# percepts per simulation:

# N looks pretty similar across noise levels and
# delta_f

p = ggplot(filter(summary,measure == "N",condition == "freqs"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(delta_f~W_m_σ,
             labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                   cols = W[m[sigma]] == .(W_m_σ))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")

save_plot(file.path(dir,"bistable_freq_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=6,base_width=2,base_height=2)

########################################
# good view of the selectivity of bistability
for (noise in unique(summary$c_σ)){
  p1 = ggplot(filter(summary,measure == "mean_ratio",condition == "freqs",
                     c_σ == noise,W_m_σ != 2),
         aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
             fill=clamp(value,-1,1))) +
    geom_raster() +
    facet_grid(delta_f~W_m_σ,
               labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                     cols = W[m[sigma]] == .(W_m_σ))) +
    scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                         breaks=c(-1,0,1),
                         labels=c("-1 (10x shorter)",
                                  " 0 (Equal Lengths)",
                                  " 1 (10x longer)")) +
    xlab(expression(paste("Adaptation ", (c[a])))) +
    ylab(expression(paste("Inhibition ", (c[m])))) +
    ggtitle(paste("Ratio of Percept Lengths (noise =",noise,")"))

  fuse_split = summary %>%
    filter(measure == "mean_ratio", condition == "freqs", c_σ == noise,
           W_m_σ != 2) %>%
    select(-pindex) %>%
    mutate(value = clamp(value,-1,1)) %>%
    spread(delta_f,value) %>%
    mutate(quality = -(abs(`3`) + abs(1-`0.5`) + abs(-1 - `12`)))

  p2 = ggplot(fuse_split,
              aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=quality)) +
    geom_raster() +
    facet_grid(~W_m_σ,labeller=
               label_bquote(cols = paste("Breadth ",(W[m[sigma]]) == .(W_m_σ)))) +
    scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                         limits=c(-3,0)) +
    xlab(expression(paste("Adaptation ", (c[a])))) +
    ylab(expression(paste("Inhibition ", (c[m])))) +
    ggtitle(expression(paste("Selectivity to ",
                             Delta[f],"= 3: -(",
                             abs(Delta[3]) - abs(1-Delta[0.5]) -
                               abs(-1 - Delta[12]),") ")))

  p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align='v')
  filename = sprintf("bistable_freq_selectivity_noise%02.2f.pdf",noise)
  save_plot(file.path(dir,filename),p,base_height=2,base_width=2,nrow=4,ncol=6)
}

########################################
# individual tracks of a single parameter
sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-42) < 6e-1,abs(c_a-32) < 6e-1,
         delta_f==3,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = c(0,cumsum(length[1:length(length)-1])))

ggplot(selected,aes(x=time,y=id,color=factor(stimulus),group=id)) +
  geom_line(size=5) +
  scale_color_brewer(palette='Set1')

#############################################################################
# plotting of object-level condition

########################################
# percepts per simulation:

# N looks pretty similar across noise levels and
# delta_f

p = ggplot(filter(summary,measure == "N",condition == "track"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,
             labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                   cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")

save_plot(file.path(dir,"bistable_track_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=6,base_width=2,base_height=2)

########################################
# selectivity of bistability

p1 = ggplot(filter(summary,measure == "mean_ratio",condition == "track"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=
             label_bquote(rows = Delta[f] == .(delta_f),
                          cols = paste("Noise ", (c[sigma] == .(c_σ))))) +
  scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                       limits=c(-1,1),
                       breaks=c(-1,0,1),
                       labels=c("-1 (10x shorter)",
                                " 0 (Equal Lengths)",
                                " 1 (10x longer)")) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle("Ratio of Percept Lengths")
p1

fuse_split = summary %>%
  filter(measure == "mean_ratio", condition == "track") %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(delta_f,value) %>%
  mutate(quality = -(abs(`3`) + abs(1-`0.5`) + abs(-1 - `12`)))

p2 = ggplot(fuse_split,
            aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=quality)) +
  geom_raster() +
  facet_grid(~c_σ,labeller=
             label_bquote(cols = paste("Noise ", (c[sigma]) == .(c_σ)))) +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-3,0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",Delta[f],"= 3: -(",
                           abs(Delta[3]) - abs(1-Delta[0.5]) -
                             abs(-1 - Delta[12]),") ")))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align="v")
save_plot(file.path(dir,"bistable_track_selectivity.pdf"),p,
          base_height=2,base_width=2,nrow=4,ncol=6)

