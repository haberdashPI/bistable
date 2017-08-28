library(ggplot2)
library(dplyr)
library(feather)
########################################
## Parameter gamma

df = read.csv('../../data/noise_level_2017-08-07_17.54.csv')

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

ggsave(paste("../../plots/noise_level_",Sys.Date(),".pdf",sep=""),
       width=11,height=8)

########################################
# Parameter sigma_p

df = read.csv('../../data/input_spread_2017-08-07_20.06.csv')

means = df %>%
  group_by(DF,PR,sigma_p) %>%
  summarize(x = mean(x))

DFs = floor(unique(means$DF) / 0.25)*0.25
PRs = floor(unique(means$PR) / 0.25)*0.25
ggplot(means,aes(x=log(DF),y=log(PR),
                 fill = cut(x,0:11/11),
                 label=round(x,2))) +
  geom_raster(interpolate=F) + geom_text(size=3) +
  facet_wrap(~paste("input spread =",sprintf("%05.2f",round(sigma_p,2)),
                    "(semitones)")) +
  scale_x_continuous(breaks=log(DFs),labels=DFs) +
  scale_y_continuous(breaks=log(PRs),labels=PRs) +
  xlab(expression(paste(Delta,"f (semitones)"))) +
  ylab("Presentation rate (Hz)") +
  scale_fill_brewer(palette='RdBu',drop=F) + theme_classic()

ggsave(paste("../../plots/input_spread_",Sys.Date(),".pdf",sep=""),
       width=11,height=8)


########################################
## phase lengths
df = read.csv('../../data/phase_lengths_2017-08-08_10.27.csv')

normalized = df %>%
  group_by(repeat.) %>%
  mutate(nlength = length/mean(length))

ggplot(subset(normalized,index > 1),aes(x=nlength)) +
  geom_histogram(bins=50) + theme_classic() + xlab('Normalized Phase Length')

ggsave(paste("../../plots/phase_lengths_",Sys.Date(),".pdf",sep=""),
       width=5,height=5)



########################################
## phase lengths
df = read.csv('../../data/int_cont_phase_lengths_2017-08-14_09.15.csv')

normalized = df %>%
  group_by(repeat.) %>%
  mutate(nlength = length/mean(length))

ggplot(subset(normalized,index > 1),aes(x=nlength)) +
  geom_histogram(bins=50) + theme_classic() + xlab('Normalized Phase Length') +
  xlim(0,8)

ggsave(paste("../../plots/int_cont_phase_lengths_",Sys.Date(),".pdf",sep=""),
       width=5,height=5)


########################################
## intermittant van noorden
df = read.csv('../../data/int_van_noorden_2017-08-11_16.06.csv')

means = df %>%
  group_by(DF,PR) %>%
  summarize(x = mean(x))

DFs = floor(unique(means$DF) / 0.25)*0.25
PRs = floor(unique(means$PR) / 0.25)*0.25
ggplot(means,aes(x=log(DF),y=log(PR),
                 fill = cut(x,0:11/11),
                 label=round(x,2))) +
  geom_raster(interpolate=F) + geom_text(size=3) +
  scale_x_continuous(breaks=log(DFs),labels=DFs) +
  scale_y_continuous(breaks=log(PRs),labels=PRs) +
  xlab(expression(paste(Delta,"f (semitones)"))) +
  ylab("Presentation rate (Hz)") +
  scale_fill_brewer(palette='RdBu',drop=F) + theme_classic()

ggsave(paste("../../plots/int_van_noorden_",Sys.Date(),".pdf",sep=""),
       width=11,height=8)

########################################
## intermittant phase lengths
df = read.csv('../../data/int_phase_lengths_2017-08-11_14.40.csv')

normalized = df %>%
  group_by(repeat.) %>%
  mutate(nlength = length/mean(length))

ggplot(subset(normalized,index > 1),aes(x=nlength)) +
  geom_histogram(bins=50) + theme_classic() + xlab('Normalized Phase Length')

ggsave(paste("../../plots/int_phase_lengths_",Sys.Date(),".pdf",sep=""),
       width=5,height=5)

##################################################
## switch input levels (in intermittent stimulus)

df = read_feather('../../data/switch_input_levels_2017-08-18_13.15.feather')
ggplot(df,aes(x=nearby_level,..density..,fill=switch)) +
  geom_histogram(alpha=0.5,bins=50)

means = NULL
for(t in seq(1e-4,0.25,length.out=100)){
  means0 = df %>%
    mutate(less_than_t = nearby_level < t) %>%
    group_by(less_than_t) %>%
    summarize(prop_switch = mean(switch))
  means0$thresh = t
  means = rbind(means,means0)
}

ggplot(means,aes(y=100*prop_switch,x=thresh,color=less_than_t)) +
  geom_line(size=1.2) +
  theme_classic(base_size=14) +
  theme(legend.position = c(0.9, 0.9),
        legend.justification='right') +
  ylab('% Switching') +
  xlab(expression(paste("Threshold (",theta,")"))) +
  scale_color_brewer("",palette='Set1',
                     labels=c(expression(paste("Recent Input",phantom(x)>=theta)),
                              expression(paste("Recent Input",phantom(x)<theta))))

ggsave(paste("../../plots/switch_input_levels_",Sys.Date(),".pdf",sep=""),
       width=6,height=4,useDingbats=F)


prop_in_gaps = df %>%
  mutate(in_gap = nearby_level < 0.01) %>%
  group_by(in_gap) %>%
  do(data.frame(rbind(smean.cl.boot(.$switch))))

ggplot(prop_in_gaps,aes(x=in_gap,y=Mean,ymin=Lower,ymax=Upper)) +
  geom_pointrange() +
  theme_classic(base_size = 14) +
  ylab('% switching') +
  xlim('During Stimulus', 'During Gap')

ks.test(subset(df,switch)$nearby_level,subset(df,!switch)$nearby_level)
