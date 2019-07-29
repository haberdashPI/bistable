
include(joinpath(@__DIR__,"interactive_setup.jl"))

datadir = joinpath(@__DIR__,"..","data","buildup")
files = readdir(datadir)

df = mapreduce(vcat,files) do file
  if !isfile(joinpath(datadir,file))
    return DataFrame()
  end
  df = DataFrame(CSV.File(joinpath(datadir,file)))
  m = match(r"^([a-z]+)_df([0-9]+)_([0-9-]+)\.csv",file)
  if isnothing(m)
    @warn("Filename $file doesn't match the expected naming convention.")
  else
    df[!,:level] .= m[1]
    df[!,:df] .= parse(Int,m[2])
    df[!,:date] .= Date(m[3])
    df
  end
end

means = by(df,[:level,:df]) do df
  buildup_mean(df,delta=0.2,length=12)
end
R"""
ggplot($means,aes(x=time,y=value,color=factor(df))) + geom_line() +
  facet_wrap(~level)
"""

obj = df[(df.level .== "object") .& (df.df .== 6),:]
img = buildup_image(df,delta=0.05,length=12)

R"""
ggplot($img,aes(x=time,y=run,fill=value)) + geom_raster()
"""
