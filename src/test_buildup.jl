using ClusterManagers
using Distributed

include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)
settings = joinpath("..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24


if endswith(gethostname(),".cluster")
    addprocs(SlurmManager(20), partition="CPU", t="04:00:00",
             cpus_per_task=1)
    @everywhere include(joinpath(@__DIR__,"..","src","setup.jl"))
end

p = copy(params[select_params(params,t_c_a=5,t_c_m=5,Δf=6),:])
p.Δf .= 6
results = @distributed (vcat) for i in 1:4000
    with_logger(NullLogger()) do
        len,val = bistable_model(p,settings,intermediate_results=true).percepts.counts
        DataFrame(length=len,response=val.+1,run=i)
    end
end

Δ = 0.05
sdf = buildup_mean(buildup,delta=Δ,length=12)

