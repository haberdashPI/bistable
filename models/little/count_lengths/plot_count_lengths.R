library(cowplot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)

dir = file.path('..','..','..','plots',paste('scale_percept_lengths_',
                                        Sys.Date(),sep='_'))
dir.create(dir,showWarnings=F)
df = read_feather(file.path('..','..','..','data','count_lengths',
                            'scale_percept_lengths_2018-05-29.feather'))
params = read_feather(file.path('..','..','..','data','count_lengths',
                                'scale_percept_params_2018-05-29.feather'))
params$pindex = 1:nrow(params)

zero.to.na = function(xs){
  ys = xs
  ys[ys == 0] = NA
  ys
}

summary = df %>% group_by(pindex,method) %>%
  summarize(N1 = sum(stimulus == 1),
	    N2 = sum(stimulus > 1),
	    N_ratio = log10((N1+1)/(N2+1)),
	    mean_length_1 = mean(log10(zero.to.na(length[stimulus == 1])),na.rm=T),
	    mean_length_2 = mean(log10(zero.to.na(length[stimulus > 1])),na.rm=T),
	    mean_ratio = mean_length_1 - mean_length_2,
	    sd_length_1 = sd(log10(zero.to.na(length[stimulus == 1])),na.rm=T),
	    sd_length_2 = sd(log10(zero.to.na(length[stimulus > 1])),na.rm=T),
	    sd_ratio = sd_length_1 - sd_length_2) %>%
  left_join(params) %>%
  gather(measure,value,N1:sd_ratio)

plot_taus = function(summary){
  ggplot(summary,aes(x = tau_m/1000,y=tau_a/1000,fill=value)) +
  facet_grid(measure~paste("noise =",tau_sigma/1000," (s)")) + geom_raster() +
  scale_x_continuous(name="Mutual Inhibition Tau (s)") +
  scale_y_continuous(name="Adaptation Tau (s)")
}

################################################################################
# By Threshold
strength_plot = function(df){
  val_breaks=c(0.2,0.5,1,2,5)
  m_breaks = c(0.1,1,10)
  a_breaks = m_breaks

  pdf = df %>%
    group_by(c_m,c_a,c_sigma,method) %>%
    summarize(value=median(value,na.rm=T))

  ggplot(pdf,aes(x = log10(c_m),y=log10(c_a),fill=value)) +
    facet_grid(method~paste("noise =",round(c_sigma,3))) + geom_raster() +
    scale_x_continuous(name="Mutual Inhibition",breaks=log10(m_breaks),labels=m_breaks) +
    scale_y_continuous(name="Adaptation",breaks=log10(a_breaks),labels=a_breaks) +
    scale_fill_distiller(name="log_10(Ratio)",palette='RdBu',
                         limits=c(-max(abs(pdf$value),na.rm=T),
                                  max(abs(pdf$value),na.rm=T))) +
    ggtitle("Mean Length Ratio (Threshold Hueristic) by Strengths")
}

p = summary %>% filter(measure == 'mean_ratio') %>% strength_plot
save_plot(file.path(dir,"mean_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
          base_height=2.5)

p = summary %>% filter(measure == 'sd_ratio') %>% strength_plot
save_plot(file.path(dir,"sd_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
                    base_height=2.5)

p = summary %>% filter(measure == 'N_ratio') %>% strength_plot
save_plot(file.path(dir,"count_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
                    base_height=2.5)

# TODO: do the same for taus
# TODO: do the same for W_m_sig and c_m (and tau_m?)
# TODO: do the same for tau_x and c_x

# TODO: though, although I was missing something before,
# this view also obscures whatever I saw in the prior plots
