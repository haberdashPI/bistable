using ClusterManagers
using Distributed
using Iterators: product

include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)
settings = joinpath(@__DIR__,"..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24
settings["bandwidth_ratio"]["window"] = 1.0
settings["bandwidth_ratio"]["delta"] = 0.25

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] = 1:size(params,1)
settings = joinpath("..","src","settings.toml")

results = []
for_results_in(joinpath(datadir,"data"),reinterpret="reinterpret") do entry
  push!(results,DataFrame(length=entry["lengths"],
                          percepts=entry["percepts"].+1, # after +1, indicates the number of streams reported, 1 or 2
                          created=entry["created"],
                          pindex=entry["pindex"])) # the parameter index (pindex = N correspondes to row N of `params`)
end
df = vcat(results...);

fields = [:f_c_a,:f_c_m,:f_c_σ,:s_c_a,:s_c_m,:s_c_σ,:t_c_a,:t_c_m,:t_c_σ]
progress = Progress(nrow(unique(params[:,fields])))
herr = human_error()
df_summary = by(params,fields) do row
    next!(progress)
    err = model_error(df,params;bound=false,(k => row[1,k] for k in fields)...)
    DataFrame(stream_error = err.stream,length_error = err.lengths,eratio = error_ratio(err,herr))
end

function level(row)
  if row.t_c_σ > 0
    @assert row.s_c_σ == 0
    @assert row.f_c_σ == 0
    "object"
  elseif row.s_c_σ > 0
    @assert row.t_c_σ == 0
    @assert row.f_c_σ == 0
    "central"
  elseif row.f_c_σ >0
    @assert row.s_c_σ == 0
    @assert row.t_c_σ == 0
    "peripheral"
  else
    error("Could not infer level")
  end
end

df_summary[!,:level] = level.(eachrow(df_summary))
thresholds = by(df_summary,:level,thresh = :eratio => x->quantile(x,0.05))

df_best = by(df_summary,:level) do group
  threshold = quantile(group.eratio,0.05)
  group[group.eratio .<= threshold,Not(:level)]
end
df_best[!,:index] = 1:size(df_best,1)

N = sum(df_best.level .== "object")
obji = df_best.index[df_best.level .== "object"]
ceni = df_best.index[df_best.level .== "central"]
peri = df_best.index[df_best.level .== "peripheral"]

combined = map(product(1:N,1:N,1:N)) do (i,j,k)
  ri = obji[i]
  rj = ceni[j]
  rk = peri[k]
  (rows = (ri,rj,rk), eratio = mean(df_best.eratio[[ri,rj,rk]]))
end |> vec
combined = sort!(combined,by=x->x.eratio)[1:N]

for comb in combined
  columns = Not([:level,:index])
  rows = df_best[collect(comb.rows),columns]
  @show rows
  newrow = (;(col => sum(rows[:,col]) for col in names(rows))...)
  push!(df_best,(level = "combined", index = size(df_best,1)+1, newrow...))
end

params = [:f_c_a,:f_c_m,:f_c_σ,:s_c_a,:s_c_m,:s_c_σ,:t_c_a,:t_c_m,:t_c_σ]

if endswith(gethostname(),".cluster")
    addprocs(SlurmManager(20), partition="CPU", t="4:00:00",
             cpus_per_task=1)
    @everywhere include(joinpath(@__DIR__,"..","src","setup.jl"))
end

writedir = joinpath(@__DIR__,"..","data","buildup")
CSV.write(joinpath(writedir,"model_params.csv"),df_best[:,Not(:index)])

N = 10^1
models = enumerate(eachrow(df_best))
deltas = [3,6,12]
@distributed (vcat) for ((mindex,model),Δ,i) in product(models,deltas,1:N)
  p = copy(AuditoryBistabilityLE.read_params,model[params])
  p.Δf .= Δ
  with_logger(NullLogger()) do
    results = bistable_model(p,settings,intermediate_results=true)
    len,val = results.percepts.counts
    DataFrame(
      length=len,
      response=val.+1,
      run=i,
      model_index = mindex,
    )
  end
end
CSV.write(joinpath(writedir,"build_results.csv"),results)
