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
                            "scale_freq_percept_lengths_2018-07-07.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "params_2018-07-03.feather"))
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

# TODO: break up by freq and scales

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

sim_len = df %>%
  filter(pindex == 1,created == first(created)) %>%
  select(length) %>% sum

summary = summary %>%
  mutate(N_per_sec_n   = norm_ratio(N / num_sims/sim_len),
         N_ratio_n     = norm_ratio(N_ratio),
         mean_ratio_n  = norm_ratio(mean_ratio),
         sd_ratio_n    = norm_ratio(sd_ratio),
         match         = pmax(W,N_ratio_n,mean_ratio_n,sd_ratio_n)) %>%
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
  filter(measure == "mean_ratio", condition == "scales") %>%
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

p1 = ggplot(filter(summary,measure == "mean_ratio",condition == "freqs"),
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
  ggtitle("Ratio of Percept Lengths")
p1


save_plot(file.path(dir,"bistable_freq_selectivity.pdf"),p1,base_aspect_ratio=1.3,
          nrow=3,ncol=6,base_width=2,base_height=2)


