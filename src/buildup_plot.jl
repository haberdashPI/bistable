
include(joinpath(@__DIR__,"interactive_setup.jl"))
curplotdir = joinpath(plotdir,string(Date(now())))
isdir(curplotdir) || mkdir(curplotdir)

datadir = joinpath(@__DIR__,"..","data","buildup","2019-08-25")
df1 = DataFrame(CSV.File(joinpath(datadir,"build_results.csv")))
df2 = DataFrame(CSV.File(joinpath(datadir,"build_results_object_only.csv")))
df = vcat(df1[in.(df1.level,Ref(["peripheral","central"])),:],df2)

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
