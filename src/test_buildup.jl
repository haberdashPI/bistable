using ClusterManagers
using Distributed

include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)
settings = joinpath(@__DIR__,"..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24


if endswith(gethostname(),".cluster")
    addprocs(SlurmManager(20), partition="CPU", t="04:00:00",
             cpus_per_task=1)
    @everywhere include(joinpath(@__DIR__,"..","src","setup.jl"))
end

# todo: this might be why we didn't get the right thing
# in our object-level simulation
params.f_c_σ .= 0.5
params.s_c_σ .= 0.5
params.t_c_σ .= 0.5

models = Dict(
  :object => copy(params[select_params(params,t_c_a=5,t_c_m=5,Δf=6),:]),
  # :central => copy(params[select_params(params,s_c_a=5,s_c_m=5,Δf=6),:]),
  # :peripheral => copy(params[select_params(params,f_c_a=15,f_c_m=130,Δf=6),:]),
)
# p = copy(models[:peripheral])
# p.f_c_σ .= 0.5
# p.s_c_a .= 5
# p.s_c_m .= 5
# p.s_c_σ .= 0.5
# p.t_c_a .= 5
# p.t_c_m .= 5
# p.t_c_σ .= 0.5
# models[:combined] = p

writedir = joinpath(@__DIR__,"..","data","buildup")
for (name,model) in models
  for Δ in [3,6,12]
    p = copy(model)
    p.Δf .= Δ
    results = @distributed (vcat) for i in 1:10^4
      with_logger(NullLogger()) do
        results = bistable_model(p,settings,intermediate_results=true)
        len,val = results.percepts.counts
        DataFrame(length=len,response=val.+1,run=i)
      end
    end
    CSV.write(joinpath(writedir,"$(name)_df$(Δ)_$(Date(now())).csv"),results)
  end
end


