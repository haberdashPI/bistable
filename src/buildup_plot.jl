
include(joinpath(@__DIR__,"interactive_setup.jl"))
curplotdir = joinpath(plotdir,string(Date(now())))
isdir(curplotdir) || mkdir(curplotdir)

datadir = joinpath(@__DIR__,"..","data","buildup","2019-08-24")
df = DataFrame(CSV.File(joinpath(datadir,"build_results.csv")))
models = DataFrame(CSV.File(joinpath(datadir,"model_params.csv")))
df[!,:level] = categorical(models.level[df.model_index])

means = by(df,[:delta_f,:level]) do df
  buildup_mean(df,delta=0.3,length=12)
end

# means = by(df,[:level,:df]) do df
#   buildup_mean(df,delta=0.3,length=12)
# end

R"""
library(dplyr)
df = $means
df$level = factor(df$level, levels=c("peripheral","central","object","combined"), ordered=T)
df = df %>% arrange(level)
ggplot(df,aes(x=time,y=value,color=factor(delta_f))) + geom_line() +
  facet_grid(~level) + xlim(0,10) + theme_classic()
ggsave($(joinpath(curplotdir,"buildup.pdf")))
"""
