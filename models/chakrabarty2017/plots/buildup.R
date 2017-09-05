require(stringr)
require(tidyr)
require(dplyr)
require(ggplot2)

## df = read.csv('../../../data/buildup.csv')
## names(df) = c('delta_1','delta_3','delta_6','delta_9','duration')
## dfr = gather(df,delta,pstreaming,delta_1:delta_9)
## dfr$delta = as.numeric(str_extract(dfr$delta,"([0-9]+)"))

## ggplot(dfr,aes(x=duration/1000,y=pstreaming,
##                group=factor(delta),color=factor(delta))) +
##   geom_line() + xlab('duration (s)') + ylab('% streaming') +
##   theme_classic(base_size = 14) +
##   scale_color_brewer(palette='Spectral',name='Delta F (st)') +
##   theme(legend.position = c(0.9,0.5))

## ggsave(paste("../../../plots/deb_buildup_",Sys.Date(),".pdf",sep=""),
##        width=6,height=4)

df = read.csv('../../../data/buildup_05-Sep-2017.csv')

ggplot(df,aes(x=time/1000,y=response,
              group=factor(delta),color=factor(delta))) +
  facet_wrap(~freq) +
  geom_line() + xlab('duration (s)') + ylab('% streaming') +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.5,0.15))
ggsave(paste("../../../plots/deb_buildup_by_freq_",Sys.Date(),".pdf",sep=""),
       width=6,height=6)

mdf = df %>%
  group_by(delta,time) %>%
  summarize(response = mean(response))

ggplot(mdf,aes(x=time/1000,y=response,
              group=factor(delta),color=factor(delta))) +
  geom_line() + xlab('duration (s)') + ylab('% streaming') +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Spectral',name='Delta F (st)') +
  theme(legend.position = c(0.15,0.8))
ggsave(paste("../../../plots/deb_buildup_",Sys.Date(),".pdf",sep=""),
       width=6,height=4)
