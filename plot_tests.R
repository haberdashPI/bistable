library(ggplot2)
library(dplyr)


########################################
## Parameter gamma

df = read.csv('noise_level_2017-08-07_17.54.csv')

means = df %>%
  group_by(DF,PR,gamma) %>%
  summarize(x = mean(x))

DFs = floor(unique(means$DF) / 0.25)*0.25
PRs = floor(unique(means$PR) / 0.25)*0.25
ggplot(means,aes(x=log(DF),y=log(PR),
                 fill = cut(x,0:11/11),
                 label=round(x,2))) +
  geom_raster(interpolate=F) + geom_text(size=3) +
  facet_wrap(~paste("gamma =",round(gamma,3))) +
  scale_x_continuous(breaks=log(DFs),labels=DFs) +
  scale_y_continuous(breaks=log(PRs),labels=PRs) +
  scale_fill_brewer(palette='RdBu') + theme_classic()
ggsave(paste("noise_level_",Sys.Date(),".pdf",sep=""),width=11,height=8)

########################################
# Parameter sigma_p

df = read.csv('input_spread_2017-08-07_18.44.csv')

means = df %>%
  group_by(DF,PR,sigma_p) %>%
  summarize(x = mean(x))

DFs = floor(unique(means$DF) / 0.25)*0.25
PRs = floor(unique(means$PR) / 0.25)*0.25
ggplot(means,aes(x=log(DF),y=log(PR),
                 fill = cut(x,0:11/11),
                 label=round(x,2))) +
  geom_raster(interpolate=F) + geom_text(size=3) +
  facet_wrap(~paste("sigma_p =",sprintf("%06.3f",round(sigma_p,3)))) +
  scale_x_continuous(breaks=log(DFs),labels=DFs) +
  scale_y_continuous(breaks=log(PRs),labels=PRs) +
  scale_fill_brewer(palette='RdBu') + theme_classic()

ggsave(paste("input_spread_",Sys.Date(),".pdf",sep=""),width=11,height=8)


########################################
## phase lengths
df = read.csv('phase_lengths_gamma0.05_2017-08-07_17.33.csv')

normalized = df %>%
  group_by(repeat.) %>%
  mutate(lengths = lengths/mean(lengths))

ggplot(normalized,aes(x=lengths)) + geom_histogram()

ggplot(normalized,aes(x=lengths,y=index)) + geom_point()
