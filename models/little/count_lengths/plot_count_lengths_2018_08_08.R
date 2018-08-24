
library(cowplot)
library(moments)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)
source("util.R")

dir = file.path("..","..","..","plots",paste("freq_percept_lengths_",
                                             Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)
df = read_feather(file.path("..","..","..","data","count_lengths",
                            "freq_percept_lengths_2018-08-14.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "freq_params_2018-08-05.feather"))
params$pindex = 1:nrow(params)

framerate = 0.25 # seconds
threshold = 2.1
bthreshold = 0.5
minlength = 0.500

times = list(component = seq(6*0.075,684*0.075,3*0.075),
             bandwidth = seq(0.5,48.0,0.25))

df_aligned = df %>% group_by(pindex,created,kind) %>%
  mutate(time = times[kind][[1]]) %>%
  do(data.frame(ratio=interpolate_times(.$ratio,.$time,times$component),
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
  filter(pindex == 1,created == first(created)) %>%
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


########################################
# percepts per simulation:

# N looks pretty similar across noise levels and
# delta_f

p = ggplot(filter(summary,measure == "N",condition == "freqs"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=value/num_sims/sim_len)) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Number of Percepts per Second")

save_plot(file.path(dir,"bistable_freq_N.png"),p,base_aspect_ratio=1.3,
          nrow=3,ncol=1,base_width=4,base_height=2)

########################################
# proporition of 1 vs. 2

p1 = ggplot(filter(summary,measure == "mean_ratio",condition == "freqs"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(Δf~.,labeller=label_bquote(rows = Delta[f] == .(Δf))) +
  scale_fill_distiller(name="Fused Percept",palette="RdBu",direction=1,
                       breaks=c(-1,0,1),
                       labels=c("-1 (10x shorter)",
                                " 0 (Equal Lengths)",
                                " 1 (10x longer)")) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(paste("Ratio of Percept Lengths"))
p1

fuse_split = summary %>%
  filter(measure == "mean_ratio", condition == "freqs") %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(Δf,value) %>%
  mutate(quality = -(abs(`3`) + abs(1-`0.5`) + abs(-1 - `12`)))

p2 = ggplot(fuse_split,
            aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=quality)) +
  geom_raster() +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1,
                       limits=c(-3,0)) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",
                           Delta[f],"= 3: -(",
                           abs(Delta[3]) - abs(1-Delta[0.5]) -
                             abs(-1 - Delta[12]),") ")))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align='v')
save_plot(file.path(dir,"bistable_freq_selectivity.png"),p,
          base_height=2,base_width=5.5,nrow=4,ncol=1)

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



