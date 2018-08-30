library(cowplot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)

dir = file.path("..","..","..","plots",paste("scale_histogram",
                                             Sys.Date(),sep="_"))
dir.create(dir,showWarnings=F)

df = read_feather(file.path("..","..","..","data","count_lengths",
                            "scale_percept_lengths_2018-05-29.feather"))
params = read_feather(file.path("..","..","..","data","count_lengths",
                                "scale_percept_params_2018-05-29.feather"))
params$pindex = 1:nrow(params)

pdf = df %>% filter(pindex == 8701)
ggplot(pdf,aes(x=length,fill=stimulus>1)) + geom_histogram()

