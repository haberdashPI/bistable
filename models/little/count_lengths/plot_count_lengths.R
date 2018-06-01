library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)

dir = file.path('..','..','plots',paste('run',Sys.Date(),sep='_'))
dir.create(dir,showWarnings=F)
df = read_feather(file.path('..','..','data','count_lengths',
			    'all_rows.feather'))
params = read_feather(file.path('..','..','data','count_lengths',
				'params.feather'))
params$pindex = 1:nrow(params)

summary = df %>% group_by(pindex,method) %>%
  summarize(N1 = sum(stimulus == 1),
	    N2 = sum(stimulus > 1),
	    N_ratio = log(N1/N2),
	    mean_length_1 = mean(log(length[stimulus == 1])),
	    mean_length_2 = mean(log(length[stimulus > 1])),
	    mean_ratio = mean_length_1 - mean_length_2,
	    sd_length_1 = sd(log(length[stimulus == 1])),
	    sd_length_2 = sd(log(length[stimulus > 1])),
	    sd_ratio = sd_length_1 - sd_length_2) %>%
  left_join(params) %>%
  gather(measure,value,N1:sd_ratio)

plot_cs = function(summary){
  m_breaks = c(0.1,1,10)
  a_breaks = m_breaks
  ggplot(summary,aes(x = log(c_m),y=log(c_a),fill=value)) +
  facet_grid(measure~paste("noise =",round(c_sigma,3))) + geom_raster() +
  scale_x_continuous(name="Mutual Inhibition",breaks=log(m_breaks),labels=m_breaks) +
  scale_y_continuous(name="Adaptation",breaks=log(a_breaks),labels=a_breaks)
}

plot_taus = function(summary){
  ggplot(summary,aes(x = tau_m/1000,y=tau_a/1000,fill=value)) +
  facet_grid(measure~paste("noise =",tau_sigma/1000," (s)")) + geom_raster() +
  scale_x_continuous(name="Mutual Inhibition Tau (s)") +
  scale_y_continuous(name="Adaptation Tau (s)")
}

val_breaks=c(0.2,0.5,1,2,5)
# plot_cs(filter(summary,measure %in% c('mean_length_1','mean_length_2'),
	       # method == 'threshold')) +
    # scale_fill_distiller(name="Time (s)",palette='Reds',
			   # breaks=log(val_breaks),labels=val_breaks) +
    # ggtitle("Mean Lengths (Threshold Hueristic) by Strengths")

plot_cs(filter(summary,measure == 'mean_ratio',method == 'threshold')) +
    scale_fill_distiller(name="Ratio",palette='RdBu') +
    ggtitle("Mean Length Ratio (Threshold Hueristic) by Strengths")
ggsave(file.path(dir,"strengths_threshold.pdf"))

val_breaks=c(0.2,0.5,1,2,5)
# plot_taus(filter(summary,measure %in% c('mean_length_1','mean_length_2'),
	       # method == 'threshold')) +
    # scale_fill_distiller(name="Time (s)",palette='Reds',
			   # breaks=log(val_breaks),labels=val_breaks) +
    # ggtitle("Mean Lengths (Threshold Hueristic) by Time Constants")

plot_taus(filter(summary,measure == 'mean_ratio',method == 'threshold')) +
    scale_fill_distiller(name="Ratio",palette='RdBu') +
    ggtitle("Mean Length Ratio (Threshold Hueristic) by Time Constants")
ggsave(file.path(dir,"taus_threshold.pdf"))

val_breaks=c(0.2,0.5,1,2,5)
# plot_cs(filter(summary,measure %in% c('mean_length_1','mean_length_2'),
	       # method == 'peaks')) +
    # scale_fill_distiller(name="Time (s)",palette='Reds',
			   # breaks=log(val_breaks),labels=val_breaks) +
    # ggtitle("Mean Lengths (Peaks Hueristic) by Strengths")

plot_cs(filter(summary,measure %in% c('mean_ratio'),method == 'peaks')) +
    scale_fill_distiller(name="Ratio",palette='RdBu') +
    ggtitle("Mean Length Ratio (Peaks Hueristic) by Strengths")
ggsave(file.path(dir,"strengths_peaks.pdf"))

val_breaks=c(0.2,0.5,1,2,5)
# plot_taus(filter(summary,measure %in% c('mean_length_1','mean_length_2'),
		 # method == 'peaks')) +
    # scale_fill_distiller(name="Time (s)",palette='Reds',
			   # breaks=log(val_breaks),labels=val_breaks) +
    # ggtitle("Mean Lengths (Peaks Hueristic) by Time Constants")

plot_taus(filter(summary,measure %in% c('mean_ratio'),method == 'peaks')) +
    scale_fill_distiller(name="Ratio",palette='RdBu') +
    ggtitle("Mean Length  Ratio (Peaks Hueristic) by Time Constants")
ggsave(file.path(dir,"taus_peaks.pdf"))

# plot_cs(filter(summary,measure %in% c('sd_length_1','sd_length_2')))
# plot_cs(filter(summary,measure %in% c('N1','N2')))

