library(cowplot)
library(moments)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)
source("util.R")

dir = file.path("..","plots",paste("individual_sim",Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)
df1 = read_feather(file.path("..","data","count_lengths","run_2018-09-10",
                            "individual_results.feather"))
df2 = read_feather(file.path("..","data","count_lengths","run_2018-09-12",
                            "individual_results.feather"))
df2$pindex = df2$pindex + max(df1$pindex)
df = rbind(df1,df2)

params1 = read_feather(file.path("..","data","count_lengths","run_2018-09-10",
                                "individual_params.feather"))
params2 = read_feather(file.path("..","data","count_lengths","run_2018-09-12",
                                "individual_extremes_params.feather"))
params1$pindex = 1:nrow(params1)
params2$pindex = 1:nrow(params2) + nrow(params1)
params = rbind(params1,params2)

framerate = 0.25 # seconds
threshold = 2.1
bthreshold = 0.5
minlength = 0.500

times = list(component = seq(0.5,47.9,0.2),
             bandwidth = seq(0.5,47.78,0.24))

max_noninf = df %>% filter(!is.infinite(ratio)) %>% select(ratio) %>% max

df_aligned = df %>% group_by(pindex,created,kind) %>%
  mutate(time = times[kind][[1]]) %>%
  do(data.frame(ratio=interpolate_times(pmin(max_noninf,.$ratio),
                                        .$time,times$component),
                time=times["component"][[1]])) %>%
  ungroup()

df_percepts = df_aligned %>%
  group_by(pindex,created) %>%
  spread(kind,ratio) %>%
  do(percept_lengths(.$component < threshold | .$bandwidth < bthreshold,
                     minlength,framerate)) %>%
  mutate(is_bound = set_bound(stimulus)) %>%
  ungroup()

sim_len = df_percepts %>%
  filter(pindex == 1,created == dplyr::first(created)) %>%
  select(length) %>% sum

# start by just plotting things as per previous approach
# to validate the results

summary = df_percepts %>% group_by(pindex) %>%
  summarize(num_sims      = length(unique(created)),
            N             = sum(!is_bound),
            W             = W(log10(length[!is_bound])),
            kurt          = kurtosis(log10(length[!is_bound])),
            skewness      = skewness(log10(length[!is_bound])),
            mean_ratio    = clean_ratio(stimulus,length,is_bound)) %>%
  left_join(params) %>%
  gather(measure,value,N:mean_ratio)

################################################################################
# freq-level

########################################
# percepts per simulation:

p = ggplot(filter(summary,measure == "N",f_c_σ > 0),
       aes(x=factor(round(f_c_a,0)),y=factor(round(f_c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")
p

save_plot(file.path(dir,"bistable_freq_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=1,base_width=4,base_height=2)

########################################
# proporition of 1 vs. 2

p1 = ggplot(filter(summary,measure == "mean_ratio",f_c_σ > 0),
       aes(x=factor(round(f_c_a,0)),y=factor(round(f_c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                       limits=c(-1,1),
                       breaks=c(-1,0,1),
                       labels=c("-1 (10x shorter)",
                                " 0 (Equal Lengths)",
                                " 1 (10x longer)")) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(paste("Ratio of Percept Lengths"))
p1

fuse_split = summary %>%
  filter(measure == "mean_ratio", f_c_σ > 0) %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(Δf,value) %>%
  mutate(quality = -(abs(`3`) + abs(1-`0.5`) + abs(-1 - `12`)))

p2 = ggplot(fuse_split,
            aes(x=factor(round(f_c_a,0)),
                y=factor(round(f_c_m,0)),fill=quality)) +
  geom_raster() +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-3,0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",
                           Delta[f],"= 3: -(",
                           abs(Delta[3]) - abs(1-Delta[0.5]) -
                             abs(-1 - Delta[12]),") ")))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align="v")
save_plot(file.path(dir,"bistable_freq_selectivity.pdf"),p,
          base_height=2,base_width=5.5,nrow=4,ncol=1)

# TODO: why the sudden change for all of the larger amounts of adaptation?

################################################################################
# scale-level

########################################
# percepts per simulation:

p = ggplot(filter(summary,measure == "N",s_c_σ > 0),
       aes(x=factor(round(s_c_a,0)),y=factor(round(s_c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")
p

save_plot(file.path(dir,"bistable_scale_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=1,base_width=4,base_height=2)

########################################
# balance of stimuli

p1 = ggplot(filter(summary,measure == "mean_ratio",s_c_σ > 0),
       aes(x=factor(round(s_c_a,0)),y=factor(round(s_c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                       limits=c(-1,1),
                       breaks=c(-1,0,1),
                       labels=c("-1 (10x shorter)",
                                " 0 (Equal Lengths)",
                                " 1 (10x longer)")) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(paste("Ratio of Percept Lengths"))
p1


fuse_split = summary %>%
  filter(measure == "mean_ratio", s_c_σ > 0) %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(Δf,value) %>%
  mutate(quality = -sqrt((`3`)^2 + (1-`0.5`)^2 + (-1 - `12`)^2))

p2 = ggplot(fuse_split,
            aes(x=factor(round(s_c_a)),
                y=factor(round(s_c_m)),fill=quality)) +
  geom_raster() +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-sqrt(3),0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",Delta[f] == 3)))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align="v")
save_plot(file.path(dir,"bistable_scale_selectivity.pdf"),p,
          base_height=2,base_width=5.5,nrow=4,ncol=1)

################################################################################
# object-level

########################################
# percepts per simulation:

p = ggplot(filter(summary,measure == "N",t_c_σ > 0),
       aes(x=factor(round(t_c_a,0)),y=factor(round(t_c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")
p

save_plot(file.path(dir,"bistable_track_N.pdf"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=1,base_width=4,base_height=2)

########################################
# balance of stimuli

p1 = ggplot(filter(summary,measure == "mean_ratio",t_c_σ > 0),
       aes(x=factor(round(t_c_a,0)),y=factor(round(t_c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                       limits=c(-1,1),
                       breaks=c(-1,0,1),
                       labels=c("-1 (10x shorter)",
                                " 0 (Equal Lengths)",
                                " 1 (10x longer)")) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(paste("Ratio of Percept Lengths"))
p1

fuse_split = summary %>%
  filter(measure == "mean_ratio", t_c_σ > 0) %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(Δf,value) %>%
  mutate(quality = -sqrt((`3`)^2 + (1-`0.5`)^2 + (-1 - `12`)^2))

p2 = ggplot(fuse_split,
            aes(x=factor(round(t_c_a)),
                y=factor(round(t_c_m)),fill=quality)) +
  geom_raster() +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-sqrt(3),0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",Delta[f] == 3)))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align="v")
save_plot(file.path(dir,"bistable_track_selectivity.pdf"),p,
          base_height=2,base_width=5.5,nrow=4,ncol=1)

# why are all of the higher adaptation levels now all unselective for 12 st?
# NOTE: it might be that I just need to use a different range of parameters now,
# e.g. small steps of adaptaiton between 0 and 5 but I still need to understand
# what is happening in the model for these cases

################################################################################
# histogram

sindex = params %>%
  filter(abs(t_c_m-100) < 1, abs(t_c_a-5) < 1,Δf==3) %>%
  select(pindex) %>% head(1) %>% first

percepts = df_percepts %>% filter(pindex == sindex,!is_bound)

ggplot(percepts,aes(x=length)) + geom_histogram(bins=12)
ggsave(file.path(dir,"track_hist_m100_a5.pdf"))

################################################################################
# old, unrevised code

# inspecting a few parameters...
sindex = params %>%
  filter(abs(c_m-42) < 6e-1,abs(c_a-6) < 6e-1,
         Δf==0.5,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first
df05 = df_aligned %>% filter(pindex == sindex,kind == "bandwidth")
ggplot(df05,aes(x=time,y=ratio,group=created)) + geom_line(alpha=0.2) +
  ylim(0,1)

sindex = params %>%
  filter(abs(c_m-42) < 6e-1,abs(c_a-6) < 6e-1,
         Δf==12,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first
df12 = df_aligned %>% filter(pindex == sindex,kind == "bandwidth")
ggplot(df12,aes(x=time,y=ratio,group=created)) + geom_line(alpha=0.2) +
  ylim(0,1)

quartz()

sindex = params %>%
  filter(abs(c_m-42) < 6e-1,abs(c_a-6) < 6e-1,
         Δf==3,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first
df3 = df_aligned %>% filter(pindex == sindex,kind == "bandwidth")
ggplot(df3,aes(x=time,y=ratio,group=created)) + geom_line(alpha=0.2) +
  ylim(0,1)

# buildup
select_params = function(params,m,a,freq){
  params %>%
    filter(abs(c_m-m) < 1,abs(c_a-a) < 1,Δf==freq,condition == "freqs") %>%
    select(pindex) %>% head(1) %>% first
}

sindices = c(select_params(params,42,0,3),select_params(params,65,6,3),
             select_params(params,100,10,3),select_params(params,65,0,3),
             select_params(params,100,6,3),select_params(params,42,6,3),
             select_params(params,65,10,3))

dfbuild = df_aligned %>% filter(pindex %in% sindices) %>%
  spread(kind,ratio) %>%
  group_by(time) %>%
  summarize(stream = mean(component < threshold | bandwidth < bthreshold))

ggplot(dfbuild,aes(x=time,y=stream)) + geom_line() +
  xlab('Time (s)') + ylab('% Streaming Responses')
ggsave(file.path(dir,"buildup_freq.png"))

dfbuild = df_aligned %>% filter(pindex %in% sindices) %>%
  spread(kind,ratio) %>%
  mutate(timewin = floor(time/0.5)) %>%
  group_by(timewin) %>%
  summarize(stream = mean(component < threshold | bandwidth < bthreshold),
            time = first(time))

ggplot(dfbuild,aes(x=time,y=stream)) + geom_line() +
  xlab('Time (s)') + ylab('% Streaming Responses') + ylim(0,1)
ggsave(file.path(dir,"buildup_freq_smooth.png"))



