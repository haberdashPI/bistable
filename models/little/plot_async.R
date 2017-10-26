library(feather)
library(dplyr)
library(ggplot2)

df = read_feather('../../data/async_noadapt_2017-09-21_13.41.feather')
df1 = read_feather('../../data/async_adapt_2017-09-21_14.08.feather')
df$adapt = 'noadapt'
df1$adapt = 'adapt'
df = rbind(df,df1)

ggplot(df,aes(x=aresp,y=bresp,color=condition)) +
  geom_point(aes(size=factor(freq))) + facet_grid(adapt~delta) +
  geom_abline(intercept=0,slope=1) + theme_classic()
ggsave(paste('../../plots/async_cont_adapt_v_none_',Sys.Date(),".pdf",sep=""))

ggplot(subset(df,adapt = 'adapt'),aes(x=aresp,y=bresp,color=condition)) +
  geom_point(aes(size=factor(freq))) + facet_grid(~delta) +
  geom_abline(intercept=0,slope=1) + theme_classic()
ggsave(paste('../../plots/async_adapt_alone_',Sys.Date(),".pdf",sep=""))
