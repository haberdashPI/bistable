using ClusterManagers
using Distributed
const product = Iterators.product

# setup environment
include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

# load in the parameters
params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)

# load and modify settings
settings = joinpath(@__DIR__,"..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24
settings["bandwidth_ratio"]["window"] = 1.0
settings["bandwidth_ratio"]["delta"] = 0.25

# read in the simulation results in `datadir`
results = []
for_results_in(joinpath(datadir,"data"),reinterpret="reinterpret") do entry
  push!(results,DataFrame(length=entry["lengths"],
                          percepts=entry["percepts"].+1, # after +1, indicates the number of streams reported, 1 or 2
                          created=entry["created"],
                          pindex=entry["pindex"])) # the parameter index (pindex = N correspondes to row N of `params`)
end
df = vcat(results...);

# load a summary of the data
fields = [:f_c_a,:f_c_m,:f_c_σ,:s_c_a,:s_c_m,:s_c_σ,:t_c_a,:t_c_m,:t_c_σ]
progress = Progress(nrow(unique(params[:,fields])))
herr = human_error()
df_summary = by(params,fields) do row
    next!(progress)
    err = model_error(df,params;bound=false,(k => row[1,k] for k in fields)...)
    DataFrame(stream_error = err.stream,length_error = err.lengths,eratio = error_ratio(err,herr))
end

# compute the corresponding level for each summary row
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

# find the top 5% performing models of each level
df_best = by(df_summary,:level) do group
  threshold = quantile(group.eratio,0.05)
  group[group.eratio .<= threshold,Not(:level)]
end
df_best[!,:index] = 1:size(df_best,1)

# use the top 5% of each level to get a putative set of ideal models
# for the combined set (may not be perfect, but should be decent).
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

# save the models used for this buildup simulation
writedir = joinpath(@__DIR__,"..","data","buildup",string(Date(now())))
isdir(writedir) || mkdir(writedir)
CSV.write(joinpath(writedir,"model_params.csv"),df_best[:,Not(:index)])

if endswith(gethostname(),".cluster")
    addprocs(SlurmManager(20), partition="CPU", t="4:00:00",
             cpus_per_task=1)
    @everywhere include(joinpath(@__DIR__,"..","src","setup.jl"))
end

# run the buildup simulation
N = 10^1
models = enumerate(eachrow(df_best))
deltas = [3,6,12]
runs = product(models,deltas,1:N)

results = @distributed (vcat) for ((mindex,model),Δ,i) in runs
  # parameter setup
  p = AuditoryBistabilityLE.read_params(params[1,:])
  for k in names(model[Not(:level)])
    p[k] = model[k]
  end
  p[:Δf] = Δ

  # run simulation
  with_logger(NullLogger()) do
    results = bistable_model(p,settings,intermediate_results=true)
    len,val = results.percepts.counts
    DataFrame(
      length=len,
      response=val.+1,
      run=i,
      model_index = mindex
    )
  end
end

CSV.write(joinpath(writedir,"build_results.csv"),results)

