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


W = function(x){
  if(length(x) > 5000){
    shapiro.test(sample(x,5000))[[1]]
  }else if(length(x) < 3){
    0
  }else if(all(x == x[1])){
    0
  }else{
    shapiro.test(x)[[1]]
  }
}

summary = df %>% group_by(pindex,method) %>%
  summarize(N1 = sum(stimulus == 1),
	    N2 = sum(stimulus > 1),
	    N_ratio = log10((N1+1)/(N2+1)),
	    mean_length_1 = mean(log10(zero.to.na(length[stimulus == 1])),na.rm=T),
	    mean_length_2 = mean(log10(zero.to.na(length[stimulus > 1])),na.rm=T),
      W = W(log10(length)),
      N = sum(stimulus > -1),
	    mean_ratio = mean_length_1 - mean_length_2,
	    sd_length_1 = sd(log10(zero.to.na(length[stimulus == 1])),na.rm=T),
	    sd_length_2 = sd(log10(zero.to.na(length[stimulus > 1])),na.rm=T),
	    sd_ratio = sd_length_1 - sd_length_2) %>%
  left_join(params) %>%
  gather(measure,value,N1:sd_ratio)

################################################################################
# By Threshold

strength_plot = function(df,uselog=T){
  val_breaks=c(0.2,0.5,1,2,5)
  m_breaks = c(0.1,1,10)
  a_breaks = m_breaks

  pdf = df %>%
    group_by(c_m,c_a,c_sigma,method) %>%
    summarize(value=median(value,na.rm=T))

  ggplot(pdf,aes(x = log10(c_m),y=log10(c_a),fill=value)) +
    facet_grid(method~paste("noise =",round(c_sigma,3))) + geom_raster() +
    scale_x_continuous(name="Mutual Inhibition",breaks=log10(m_breaks),
                       labels=m_breaks) +
    scale_y_continuous(name="Adaptation",breaks=log10(a_breaks),labels=a_breaks) +
    if(uselog){
      scale_fill_distiller(name="log_10(Ratio)",palette='RdBu',
                           limits=c(-max(abs(pdf$value),na.rm=T),
                                    max(abs(pdf$value),na.rm=T)))
    }
}

p = summary %>% filter(measure == 'mean_ratio') %>% strength_plot
p = p + ggtitle("Strengths - Log Ratio of Means")
save_plot(file.path(dir,"mean_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
          base_height=2.5)

p = summary %>% filter(measure == 'W') %>% strength_plot(uselog=F)
p = p + ggtitle("Strengths - Shapiro-Wilk Statistic")
save_plot(file.path(dir,"shapiro_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
          base_height=2.5)

# hold on a second, why are the aboslute N's between 10 and 20?
# is this right?? if so something about the saved data is wrong
# if that's the case, then re-run for one parameter test
# and verify that new code records more than 10
p = summary %>% filter(measure == 'N') %>% strength_plot(uselog=F)
p = p + ggtitle("Strengths - Counts")
save_plot(file.path(dir,"absolute_N_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
          base_height=2.5)

p = summary %>% filter(measure == 'sd_ratio') %>% strength_plot
p = p + ggtitle("Strengths - Log Ratio of SDs")
save_plot(file.path(dir,"sd_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
                    base_height=2.5)

p = summary %>% filter(measure == 'N_ratio') %>% strength_plot
p = p + ggtitle("Strengths - Log Ratio of Counts")
save_plot(file.path(dir,"count_strengths.pdf"),p,nrow=2,ncol=5,base_width=2,
                    base_height=2.5)


p = summary %>% filter(measure == 'N_ratio',method == 'peaks') %>% strength_plot
p = p + ggtitle("Strengths - Log Ratio of Counts")
save_plot(file.path(dir,"count_strengths_peaks.pdf"),p,nrow=1,ncol=5,base_width=2,
                    base_height=2.5)

tau_plot = function(df){
  pdf = df %>%
    group_by(tau_m,tau_a,tau_sigma,method) %>%
    summarize(value=median(value,na.rm=T))

  ggplot(pdf,aes(x = tau_m/1000,y=tau_a/1000,fill=value)) +
    facet_grid(method~paste("noise tau (s) =",round(tau_sigma/1000))) +
    geom_raster() +
    scale_x_continuous(name="Mutual Inhibition Tau (s)") +
    scale_y_continuous(name="Adaptation Tau (s)") +
    scale_fill_distiller(name="log_10(Ratio)",palette='RdBu',
                         limits=c(-max(abs(pdf$value),na.rm=T),
                                  max(abs(pdf$value),na.rm=T)))
}


p = summary %>% filter(measure == 'mean_ratio') %>% tau_plot
p = p + ggtitle("Time Constants - Log Ratio of Means")
save_plot(file.path(dir,"mean_time_constants.pdf"),p,nrow=2,ncol=2,base_width=3,
          base_height=2.5)


p = summary %>% filter(measure == 'sd_ratio') %>% tau_plot
p = p + ggtitle("Time Constants - Log Ratio of SD")
save_plot(file.path(dir,"sd_time_constants.pdf"),p,nrow=2,ncol=2,base_width=3,
          base_height=2.5)

p = summary %>% filter(measure == 'N_ratio') %>% tau_plot
p = p + ggtitle("Time Constants - Log Ratio of Counts")
save_plot(file.path(dir,"count_time_constants.pdf"),p,nrow=2,ncol=2,base_width=3,
          base_height=2.5)

p = summary %>% filter(measure == 'N_ratio',method == 'peaks') %>% tau_plot
p = p + ggtitle("Time Constants - Log Ratio of Counts")
save_plot(file.path(dir,"count_time_constants_peaks.pdf"),p,nrow=1,ncol=2,
          base_width=3,base_height=2.5)

width_plot = function(df){
  pdf = df %>%
    group_by(W_m_sig,c_m,method) %>%
    summarize(value=median(value,na.rm=T))

  m_breaks = c(0.1,1,10)
  w_breaks = c(0.1,1,10,100)
  ggplot(pdf,aes(x = log10(W_m_sig),y = log10(c_m),fill=value)) +
    facet_wrap(~method) +
    geom_raster() +
    scale_x_continuous(name="Mutual Inhibition Width",breaks=log10(w_breaks),
                       labels=w_breaks) +
    scale_y_continuous(name="Mutual Inhibition",breaks=log10(m_breaks),
                       labels=m_breaks) +
    scale_fill_distiller(name="log_10(Ratio)",palette='RdBu',
                         limits=c(-max(abs(pdf$value),na.rm=T),
                                  max(abs(pdf$value),na.rm=T)))
}

p = summary %>% filter(measure == 'mean_ratio') %>% width_plot
p = p + ggtitle("Inhibition Width - Log Ratio of Means")
save_plot(file.path(dir,"inhibit_width_mean.pdf"),p,base_aspect_ratio=1.3,
          ncol=2)

p = summary %>% filter(measure == 'sd_ratio') %>% width_plot
p = p + ggtitle("Inhibition Width - Log Ratio of SDs")
save_plot(file.path(dir,"inhibit_width_sd.pdf"),p,base_aspect_ratio=1.3,
          ncol=2)

p = summary %>% filter(measure == 'N_ratio') %>% width_plot
p = p + ggtitle("Inhibition Width - Log Ratio of Counts")
save_plot(file.path(dir,"inhibit_width_count.pdf"),p,base_aspect_ratio=1.3,
          ncol=2)

# TODO: though, although I was missing something before,
# this view also obscures whatever I saw in the prior plots
