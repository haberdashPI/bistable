
# using ClusterManagers
# using Distributed
const product = Iterators.product

include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)
settings = joinpath(@__DIR__,"..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24

settings["bandwidth_ratio"]["window"] = 1.0
settings["bandwidth_ratio"]["delta"] = 0.25

# TODO: get the parameters for the top 10 object models
function only(x)
  @assert length(x) == 1
  first(x)
end
model_ranks = 1:10
models = begin
    select = DataFrame(CSV.read(joinpath(datadir,"model_rankings.csv")) )
    sort!(select,:eratio)
    select = select[(select.t_c_m .> 0) .| (select.t_c_a .> 0),:]

    indices = [
      only(select_params(params,
        t_c_m=select[i,:t_c_m],
        t_c_a=select[i,:t_c_a],Δf=6))
      for i in model_ranks
    ]
    params[indices,:]
end

writedir = joinpath(@__DIR__,"..","data","buildup_smooth",string(Date(now())))
isdir(writedir) || mkdir(writedir)

N = 2
runs = collect(product(eachrow(models),1:N))

#@distributed (vcat) for ((name,model),Δ,i) in runs
results = mapreduce(vcat,runs) do (model,i)
  with_logger(NullLogger()) do
    results = bistable_model(model,settings,intermediate_results=true)
    len,val = results.percepts.counts
    DataFrame(
      length=len,
      response=val.+1,
      run=i,
    )
  end
end

CSV.write(joinpath(writedir,"build_object_find_smooth.csv"),results)
