
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

df = read.csv("exp1A.csv")

df = df %>%
  group_by(sid) %>%
  mutate(trial = cumsum(lag(sample,default=first(sample)) > sample))


# compute percept lens
lens = df %>%
  group_by(sid,trial) %>%
  mutate(percept = cumsum(lag(response,
                              default=first(response)) != response)) %>%
  filter(same_freq,st > 1) %>%
  group_by(sid,trial,percept) %>%
  filter(time < 6.72) %>%
  summarize(context_A = first(context_A),context_B = first(context_B),
            test_A = first(test_A), test_B = first(test_B), st = first(st),
            len = last(time) - first(time),stimulus = first(response))

ggplot(lens,aes(x=len)) + geom_histogram(bins=15)

ggplot(lens,aes(x=len,fill=factor(stimulus))) +
  geom_histogram(bins=15)

df = read.csv("exp1A.csv")
stream_prop = df %>%
  filter(st > 1,time > 3,time < 6.72,same_freq) %>%
  group_by(sid,st) %>%
  summarize(response = mean(response))

ggplot(stream_prop,aes(x=factor(st),y=response)) +
  stat_summary(geom='pointrange',width=0.1,fun.data=mean_cl_boot,
               fun.args=list(conf.int=0.95)) +
  stat_summary(geom='line',fun.data=mean_cl_boot,group=1) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic()

