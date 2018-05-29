library(tidyr)
library(dplyr)
library(ggplot2)
library(feather)

df = read_feather(file.path('..','..','data','count_lengths',
			    'all_rows.feather'))
params = read_feather(file.path('..','..','data','count_lengths',
				'params.feather'))
params$pindex = 1:nrow(params)

summary = df %>% group_by(pindex) %>%
	summarize(N1 = sum(stimulus == 1),
		  N2 = sum(stimulus > 1),
		  mean_length_1 = mean(log(length[stimulus == 1])),
		  mean_length_2 = mean(log(length[stimulus > 1])),
		  sd_length_1 = sd(log(length[stimulus == 1])),
		  sd_length_2 = sd(log(length[stimulus > 1]))) %>%
	left_join(params)





