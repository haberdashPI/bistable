require(stringr)
require(tidyr)
require(dplyr)
require(ggplot2)

df = read.csv('../../../data/alter_v_sync.csv')
names(df) = c('alternating','synchronous','delta')
dfr = gather(df,condition,dprime,alternating,synchronous)

ggplot(dfr,aes(x=delta,y=dprime,group=condition,color=condition)) +
  geom_line() + xlab('delta (st)') + ylab("d'") +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette='Set1',name='Condition') +
  theme(legend.position = c(0.15,0.9))

ggsave(paste("../../../plots/deb_alt_sync_",Sys.Date(),".pdf",sep=""),
       width=6,height=4)
