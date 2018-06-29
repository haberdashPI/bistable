library(cowplot)
library(moments)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)

dir = file.path("..","..","..","plots",paste("scale_percept_lengths_",
                                             Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)
df = read_feather(file.path("..","..","..","data","count_lengths",
                            "scale_percept_lengths_2018-06-29.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "scale_percept_params_2018-06-28.feather"))
params$pindex = 1:nrow(params)

W = function(x){
  if (sum(!is.na(x)) > 5000){
    shapiro.test(sample(x,5000))[[1]]
  }else if (sum(!is.na(x)) < 3){
    0
  }else if (all(is.na(x)) || all(x == first(x[!is.na(x)]),na.rm=T)){
    0
  }else{
    shapiro.test(x)[[1]]
  }
}

trimby = function(xs,lengths){
  xs[cumsum(lengths) < 2] = NA
  xs[1] = NA
  xs[length(xs)] = NA
  xs
}

cleaned_df = df %>% group_by(pindex,created) %>%
  mutate(length = trimby(length,length),
         stimulus = trimby(stimulus,length))

summary = cleaned_df %>% group_by(pindex) %>%
  summarize(num_sims      = length(unique(created)),

            N             = sum(!is.na(stimulus)),
            N1            = sum(stimulus == 0,na.rm=T),
            N2            = sum(stimulus == 1,na.rm=T),
            N_ratio       = log10( (N1+1) / (N2+1) ),

            W             = W(log10(length)),
            kurt          = kurtosis(log10(length),na.rm=T),
            skewness      = skewness(log10(length),na.rm=T),

            mean_length   = mean(log10(length),na.rm=T),
            mean_length_1 = mean(log10(length[stimulus == 0]),na.rm=T),
            mean_length_2 = mean(log10(length[stimulus == 1]),na.rm=T),
            mean_ratio    = mean_length_1 - mean_length_2,

            sd_length     = sd(log10(length),na.rm=T),
            sd_length_1   = sd(log10(length[stimulus == 0]),na.rm=T),
            sd_length_2   = sd(log10(length[stimulus == 1]),na.rm=T),
            sd_ratio      = sd_length_1 - sd_length_2,
            sd_ratio_n    = abs(sd_ratio / quantile(sd_ratio,0.95,na.rm=T)))

norm_ratio = function(x){
  span = quantile(x[!is.infinite(x)],c(0.05,0.95),na.rm=T)
  (x / pmax(1e-8,span[2] - span[1],na.rm=T))
}

clamp = function(x,lower,upper){
  pmin(upper,pmax(lower,x))
}

summary = summary %>%
  mutate(N_ratio_n     = norm_ratio(N_ratio),
         mean_ratio_n  = norm_ratio(mean_ratio),
         sd_ratio_n    = norm_ratio(sd_ratio),
         match         = pmax(W,N_ratio_n,mean_ratio_n,sd_ratio_n)) %>%
  left_join(params) %>%
  gather(measure,value,N:match)

################################################################################

# N looks pretty similar across noise levels and
# delta_f
ggplot(filter(summary,measure == "N"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(c[a])) + ylab(expression(c[m])) +
  ggtitle("Total Number of Percepts")

N_agg = summary %>%
  filter(measure == "N") %>%
  group_by(c_m,c_a) %>%
  summarize(N_per_sim = mean(value/num_sims))

sim_len = df %>%
  filter(pindex == 1,created == first(created)) %>%
  select(length) %>% sum

p = ggplot(N_agg,
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=N_per_sim)) +
  geom_raster() +
  scale_fill_distiller(name="N",palette="Spectral") +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle("Total Number of Percepts per Simulation",
          subtitle=paste("Simulation Length =",round(sim_len),"s"))
save_plot(file.path(dir,"bistable_scale_N.pdf"),p,base_aspect_ratio=1.3)

################################################################################
# good view of the selectivity of bistability

p1 = ggplot(filter(summary,measure == "mean_ratio"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=clamp(value,-1,1))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=
             label_bquote(rows = Delta[f] == .(delta_f),
                          cols = paste("Noise ",(c[sigma] == .(c_σ))))) +
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
  filter(measure == "mean_ratio") %>%
  select(-pindex) %>%
  mutate(value = clamp(value,-1,1)) %>%
  spread(delta_f,value) %>%
  mutate(quality = -(abs(`6`) + abs(1-`2`) + abs(-1 - `12`)))

p2 = ggplot(fuse_split,
            aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=quality)) +
  geom_raster() +
  facet_grid(~c_σ,labeller=label_bquote(cols = paste("Noise ",(c[sigma]) == .(c_σ)))) +
  scale_fill_distiller(name="Selectivity",palette="Greens",direction=1) +
  xlab(expression(paste("Adaptation ", (c[a])))) +
  ylab(expression(paste("Inhibition ", (c[m])))) +
  ggtitle(expression(paste("Selectivity to ",Delta[f],"= 6: -(",
                           abs(Delta[6]) - abs(1-Delta[2]) -
                             abs(-1 - Delta[12]),") ")))

p = plot_grid(p1,p2,nrow=2,rel_heights=c(0.65,0.35),align='v')
save_plot(file.path(dir,"bistable_scales_selectivity.pdf"),p,
          base_height=2,base_width=3,nrow=4,ncol=6)

################################################################################
# how log-normal are the data?

# TODO: summarize over different delta_t's (by max)
# TODO: create an aggregate measure

ggplot(filter(summary,measure == "W"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Wilks-Shapiro",palette="Greens",direction=1) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = summary %>% filter(measure == "kurt") %>%
  select(value) %>% abs %>% max(na.rm=T)
ggplot(filter(summary,measure == "kurt"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Kurtosis",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = summary %>% filter(measure == "skewness") %>%
  select(value) %>% abs %>% max(na.rm=T)
ggplot(filter(summary,measure == "skewness"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=value)) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Skewness",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

################################################################################
# let's plot a histogram of some parametesr
sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-27) < 6e-1,abs(c_a-6) < 6e-1,delta_f==6) %>%
  select(pindex) %>% first
selected = df %>% filter(pindex == sindex)
ggplot(selected,aes(x=length)) + geom_histogram(bins=10)
ggplot(selected,aes(x=length,fill=factor(stimulus))) + geom_histogram(bins=10) +
  scale_fill_brewer(palette='Set1')

sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-42) < 6e-1,abs(c_a-6) < 6e-1,delta_f==6) %>%
  select(pindex) %>% first
selected = df %>% filter(pindex == sindex)
ggplot(selected,aes(x=length)) + geom_histogram(bins=10)
ggplot(selected,aes(x=length,fill=factor(stimulus))) + geom_histogram(bins=10) +
  scale_fill_brewer(palette='Set1')

sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-42) < 6e-1,abs(c_a-6) < 6e-1,delta_f==12) %>%
  select(pindex) %>% first
selected = df %>% filter(pindex == sindex)
ggplot(selected,aes(x=length)) + geom_histogram(bins=10)

ggplot(selected,aes(x=length,fill=factor(stimulus))) + geom_histogram(bins=10) +
  scale_fill_brewer(palette='Set1')

sindex = params %>%
  filter(abs(c_σ-0) < 1e-1,abs(c_m-42) < 6e-1,abs(c_a-6) < 6e-1,delta_f==6) %>%
  select(pindex) %>% first
selected = df %>% filter(pindex == sindex)
ggplot(selected,aes(x=length)) + geom_histogram(bins=10)

ggplot(selected,aes(x=length,fill=factor(stimulus))) + geom_histogram(bins=10) +
  scale_fill_brewer(palette='Set1')
################################################################################

# TASKS:

# 1. email mounya graphs
# 2. test/run frequency-level bistability
# 3. test/run multi-track bistability
# 4. test/run more specific scale bistability?

# Questions:

# 3. is there a good measure (something about SD?) to
# exclude the zero noise condition?
