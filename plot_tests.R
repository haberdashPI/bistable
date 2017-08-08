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
  facet_wrap(~paste("noise SD =",round(gamma,3))) +
  scale_x_continuous(breaks=log(DFs),labels=DFs) +
  scale_y_continuous(breaks=log(PRs),labels=PRs) +
  scale_fill_brewer(palette='RdBu') + theme_classic() +
  xlab(expression(paste(Delta,"f (semitones)"))) +
  ylab("Presentation rate (Hz)")

ggsave(paste("noise_level_",Sys.Date(),".pdf",sep=""),width=11,height=8)

########################################
# Parameter sigma_p

df = read.csv('input_spread_2017-08-07_20.06.csv')

means = df %>%
  group_by(DF,PR,sigma_p) %>%
  summarize(x = mean(x))

DFs = floor(unique(means$DF) / 0.25)*0.25
PRs = floor(unique(means$PR) / 0.25)*0.25
ggplot(means,aes(x=log(DF),y=log(PR),
                 fill = cut(x,0:11/11),
                 label=round(x,2))) +
  geom_raster(interpolate=F) + geom_text(size=3) +
  facet_wrap(~paste("input spread =",sprintf("%05.2f",round(sigma_p,2)),"(semitones)")) +
  scale_x_continuous(breaks=log(DFs),labels=DFs) +
  scale_y_continuous(breaks=log(PRs),labels=PRs) +
  xlab(expression(paste(Delta,"f (semitones)"))) +
  ylab("Presentation rate (Hz)") +
  scale_fill_brewer(palette='RdBu',drop=F) + theme_classic()

ggsave(paste("input_spread_",Sys.Date(),".pdf",sep=""),width=11,height=8)


########################################
## phase lengths
df = read.csv('phase_lengths_2017-08-08_10.27.csv')

normalized = df %>%
  group_by(repeat.) %>%
  mutate(nlengths = lengths/mean(lengths))

ggplot(subset(normalized,index > 1),aes(x=nlengths)) +
  geom_histogram(bins=50) + theme_classic() + xlab('Normalized Phase Length (s)')

ggsave(paste("phase_lengths_",Sys.Date(),".pdf",sep=""),width=5,height=5)
