library(cowplot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)

dir = file.path("..","..","..","plots",paste("scale_percept_lengths_",
                                             Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)
df = read_feather(file.path("..","..","..","data","count_lengths",
                            "scale_percept_lengths_2018-06-28.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "scale_percept_params_2018-06-28.feather"))
params$pindex = 1:nrow(params)

zero_to_na = function(xs){
  ys = xs
  ys[ys == 0] = NA
  ys
}


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

  mutate(length = trimby(length,length),
         stimulus = trimby(stimulus,length))

summary = cleaned_df %>% group_by(pindex) %>%
  summarize(N             = sum(!is.na(stimulus)),
            N1            = sum(stimulus == 0,na.rm=T),
            N2            = sum(stimulus == 1,na.rm=T),
            N_ratio       = log10( (N1+1) / (N2+1) ),

            W             = W(log10(length)),

            mean_length   = mean(log10(length),na.rm=T),
            mean_length_1 = mean(log10(length[stimulus == 0]),na.rm=T),
            mean_length_2 = mean(log10(length[stimulus == 1]),na.rm=T),
            mean_ratio    = mean_length_1 - mean_length_2,

            sd_length_1   = sd(log10(length[stimulus == 0]),na.rm=T),
            sd_length_2   = sd(log10(length[stimulus == 1]),na.rm=T),
            sd_ratio      = sd_length_1 - sd_length_2,

            match = sqrt(((1-W))^2+(mean_ratio)^2+(sd_ratio)^2+(N_ratio)^2)

            ) %>%
  left_join(params) %>%
  gather(measure,value,N:match)

################################################################################
# By Threshold

zero_na = function(x){
  x[x == 0] = NA
  x
}

ggplot(filter(summary,measure == "match"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=zero_na(value)))+
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Match",palette="Spectral",direction=-1,
                       values=c(0,0.3,0.6,0.7,0.8,0.9,1)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

ggplot(filter(summary,measure == "W"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=zero_na(value)))+
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="W",palette="Spectral",direction=-1,
                       values=c(0,0.3,0.6,0.7,0.8,0.9,1)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

ggplot(filter(summary,measure == "N"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=zero_na(value))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="N",palette="Spectral",direction=-1) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = max(abs(filter(summary,measure=="N_ratio",c_a > 0,c_m > 0)$value),
             na.rm=T)
ggplot(filter(summary,measure == "N_ratio"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=zero_na(value))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="N ratio",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

bounds = range(exp(filter(summary,measure=="mean_length",
                          c_a > 0,c_m > 0)$value),na.rm=T)
ggplot(filter(summary,measure == "mean_length"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),
           fill=exp(zero_na(value)))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Mean Length (s)",palette="Spectral",direction=-1,
                       limits=bounds) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = max(abs(filter(summary,measure=="mean_ratio",c_a > 0,c_m > 0)$value),
             na.rm=T)
ggplot(filter(summary,measure == "mean_ratio"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=zero_na(value))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="Mean Ratio",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))

maxabs = max(abs(filter(summary,measure=="sd_ratio",c_a > 0,c_m > 0)$value),
             na.rm=T)
ggplot(filter(summary,measure == "sd_ratio"),
       aes(x=factor(round(c_a,0)),y=factor(round(c_m,0)),fill=zero_na(value))) +
  geom_raster() +
  facet_grid(delta_f~c_σ,labeller=label_bquote(rows = Delta[f] == .(delta_f),
                                               cols = c[sigma] == .(c_σ))) +
  scale_fill_distiller(name="SD Ratio",palette="RdBu",direction=1,
                       limits=c(-maxabs,maxabs)) +
  xlab(expression(c[a])) + ylab(expression(c[m]))
