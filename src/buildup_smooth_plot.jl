
include(joinpath(@__DIR__,"interactive_setup.jl"))
curplotdir = joinpath(plotdir,string(Date(now())))
isdir(curplotdir) || mkdir(curplotdir)

datadir = joinpath(@__DIR__,"..","data","buildup_smooth","2020-01-23")
files = readdir(datadir)
df = DataFrame(CSV.read(joinpath(datadir,files[1])))

means = buildup_mean(df,delta=0.3,length=12)

R"""
library(dplyr)
df = $means
ggplot(df,aes(x=time,y=value)) + geom_line() +
  xlim(0,10) + theme_classic()
ggsave($(joinpath(curplotdir,"buildup.pdf")))
"""
