library(feather)
library(dplyr)
library(ggplot2)

df = read_feather('../../data/buildup_cont_2017-09-05_11.46.feather')

mdf = df %>%
  group_by(delta,time) %>%
  summarize(response = mean(response))

ggplot(mdf,aes(x=time,y=response,color=factor(delta))) + geom_line() +
  xlab('duration (s)') + ylab('% streaming') +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.9,0.5))
ggsave(paste('../../plots/buildup_cont_',Sys.Date(),".pdf",sep=""))

ggplot(df,aes(x=time,y=response,color=factor(delta))) + geom_line() +
  facet_wrap(~freq) +
  xlab('duration (s)') + ylab('% streaming') +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.5,0.15))
ggsave(paste('../../plots/buildup_cont_by_freq_',Sys.Date(),".pdf",sep=""))

df = read_feather('../../data/buildup_cont_L1_2017-09-05_12.14.feather')


mdf = df %>%
  group_by(delta,time) %>%
  summarize(response = mean(response))

ggplot(mdf,aes(x=time,y=response,color=factor(delta))) + geom_line() +
  xlab('duration (s)') + ylab('% streaming') +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.9,0.5))
ggsave(paste('../../plots/buildup_cont_L1_',Sys.Date(),".pdf",sep=""))

ggplot(df,aes(x=time,y=response,color=factor(delta))) + geom_line() +
  facet_wrap(~freq) +
  xlab('duration (s)') + ylab('% streaming') +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.5,0.15))
ggsave(paste('../../plots/buildup_cont_by_freq_L1_',Sys.Date(),".pdf",sep=""))


df = read_feather('../../data/buildup_cont_adapt_2017-09-21_14.16.feather')

mdf = df %>%
  group_by(delta,time) %>%
  summarize(response = mean(response))

df = read_feather('../../data/buildup_cont_noadapt_2017-09-21_13.26.feather')

mdf2 = df %>%
  group_by(delta,time) %>%
  summarize(response = mean(response))

mdf$condition = 'adapt'
mdf2$condition = 'noadapt'
mdf = rbind(mdf,mdf2)

ggplot(mdf,aes(x=time,y=response,color=factor(delta))) + geom_line() +
  xlab('duration (s)') + ylab('% streaming') +
  facet_wrap(~condition) +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.9,0.5))

ggsave(paste('../../plots/buildup_cont_adapt_v_none_',Sys.Date(),".pdf",sep=""))
