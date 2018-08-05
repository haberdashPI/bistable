library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(feather)
library(RColorBrewer)
source("util.R")


dir = file.path("..","..","..","plots",paste("freq_ind",Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)

df = read_feather(file.path("..","..","..","data","count_lengths",
                            "scale_percept_lengths_2018-07-21.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "params_2018-07-16.feather"))
params$pindex = 1:nrow(params)

########################################
# what's going on with the fusing 12 st percepts?

sindex = params %>%
  filter(abs(c_σ-0.4) < 1e-1,abs(c_m-27) < 6e-1,abs(c_a-10) < 6e-1,
         W_m_σ == 5.6,delta_f==12,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = first(.$stimulus),id = first(.$id)))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(stimulus),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"))

########################################
# what's going on with the same parameters at 3 st?
sindex = params %>%
  filter(abs(c_σ-0.4) < 1e-1,abs(c_m-27) < 6e-1,abs(c_a-10) < 6e-1,
         W_m_σ == 5.6,delta_f==3,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = first(.$stimulus),id = first(.$id)))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(stimulus),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"))

########################################
# what's going on with the splitting 12 st percepts

sindex = params %>%
  filter(abs(c_σ-0.4) < 1e-1,abs(c_m-0) < 6e-1,abs(c_a-56) < 6e-1,
         delta_f==12,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = first(.$stimulus),id = first(.$id)))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(stimulus),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"))

## basic account
# with low inhibition and high adaptation
# both tones are present, and so a splitting percept is reportd
#
# with high inhibition and low adaptaiton, one tone wins
# and so a fused percept is reported.

########################################
# plot of a good 3st pattern

sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-65) < 6e-1,abs(c_a-6) < 6e-1,
         W_m_σ == 5.6,delta_f==3,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = 0,id = first(.$id),
             created = first(.$created))) %>%
  group_by(created) %>%
  arrange(id,time) %>%
  mutate(next_stim = lead(stimulus,default=0))

# select_test = df %>% filter(pindex == sindex) %>%
#   mutate(.,id = group_indices(.,created)) %>%
#   group_by(created) %>%
#   mutate(is_bound = set_bound(length))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(next_stim),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"),
                     name="Response") +
  ylab("Simulation #") + xlab("Time (s)")
ggsave(file.path(dir,"freq_selective_responses_3st.pdf"))

########################################
# plot of the same result for 12st (currently using the wrong metric)

sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-65) < 6e-1,abs(c_a-6) < 6e-1,
         W_m_σ == 5.6,delta_f==12,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = 0,id = first(.$id),
             created = first(.$created))) %>%
  group_by(created) %>%
  arrange(id,time) %>%
  mutate(next_stim = lead(stimulus,default=0))

# select_test = df %>% filter(pindex == sindex) %>%
#   mutate(.,id = group_indices(.,created)) %>%
#   group_by(created) %>%
#   mutate(is_bound = set_bound(length))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(next_stim),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"),
                     name="Response") +
  ylab("Simulation #") + xlab("Time (s)")
ggsave(file.path(dir,"freq_selective_responses_12st.pdf"))

########################################
# plot of the same result for 0.5st (currently using the wrong metric)

sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-65) < 6e-1,abs(c_a-6) < 6e-1,
         W_m_σ == 5.6,delta_f==0.5,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = 0,id = first(.$id),
             created = first(.$created))) %>%
  group_by(created) %>%
  arrange(id,time) %>%
  mutate(next_stim = lead(stimulus,default=0))

# select_test = df %>% filter(pindex == sindex) %>%
#   mutate(.,id = group_indices(.,created)) %>%
#   group_by(created) %>%
#   mutate(is_bound = set_bound(length))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(next_stim),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"),
                     name="Response") +
  ylab("Simulation #") + xlab("Time (s)")

ggsave(file.path(dir,"freq_selective_responses_0.5st.pdf"))



## basic acount
# when the resposnes are closer there is less competition becuase of the
# wide-spread of activation required which allows both units to respond,
# ultimately leading to a generation of two percepts

########################################
# highly selective frequency parameters

sindex = params %>%
  filter(abs(c_σ-0.2) < 1e-1,abs(c_m-18) < 6e-1,abs(c_a-56) < 6e-1,
         W_m_σ == 19.4,delta_f==3,condition == "freqs") %>%
  select(pindex) %>% head(1) %>% first

selected = df %>% filter(pindex == sindex) %>%
  mutate(.,id = group_indices(.,created)) %>%
  group_by(created) %>%
  mutate(time = cumsum(length)) %>%
  select(time,stimulus,id) %>%
  do(add_row(.,time = 0,stimulus = 0,id = first(.$id),
             created = first(.$created))) %>%
  group_by(created) %>%
  arrange(id,time) %>%
  mutate(next_stim = lead(stimulus,default=0))

colors = brewer.pal(3,"Set1")
ggplot(selected,aes(x=time,y=id,color=factor(next_stim),group=id)) +
  geom_line(size=3) +
  scale_color_manual(values=colors[c(2,1)],labels=c("1-Stream","2-Streams"))

ggsave(file.path(dir,"freq_runs_6st_m18_a56_sig_19.4.pdf"))

