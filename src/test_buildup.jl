using ClusterManagers
using Distributed

include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)
settings = joinpath(@__DIR__,"..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24
settings["bandwidth_ratio"]["window"] = 1.2
settings["bandwidth_ratio"]["delta"] = 0.3

if endswith(gethostname(),".cluster")
    addprocs(SlurmManager(20), partition="CPU", t="24:00:00",
             cpus_per_task=1)
    @everywhere include(joinpath(@__DIR__,"..","src","setup.jl"))
end

models = Dict(
  :object => begin
    p = copy(params[select_params(params,t_c_a=5,t_c_m=5,Δf=6),:])
    p.t_c_σ .= 0.5
    p
  end,
  :central => begin
    p = copy(params[select_params(params,s_c_a=5,s_c_m=5,Δf=6),:])
    p.s_c_σ .= 0.5
    p
  end,
  :peripheral => begin
    p = copy(params[select_params(params,f_c_a=15,f_c_m=130,Δf=6),:])
    p.f_c_σ .= 0.5
    p
  end,
  :combined => begin
    p = copy(params[select_params(params,f_c_a=15,f_c_m=130,Δf=6),:])
    p.f_c_σ .= 0.5
    p.s_c_a .= 5
    p.s_c_m .= 5
    p.s_c_σ .= 0.5
    p.t_c_a .= 5
    p.t_c_m .= 5
    p.t_c_σ .= 0.5
    p
  end
)

# N = 10^4
N = 1000
writedir = joinpath(@__DIR__,"..","data","buildup")
for (name,model) in models
  for Δ in [3,6,12]
    p = copy(model)
    p.Δf .= Δ
    results = mapreduce(vcat,1:N) do i #@distributed (vcat) for i in 1:N
      with_logger(NullLogger()) do
        results = bistable_model(p,settings,intermediate_results=true)
        len,val = results.percepts.counts
        DataFrame(length=len,response=val.+1,run=i)
      end
    end
    CSV.write(joinpath(writedir,"$(name)_df$(Δ)_$(Date(now())).csv"),results)
  end
end


