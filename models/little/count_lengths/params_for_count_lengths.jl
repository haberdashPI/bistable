using FileIO
using DataFrames
using Feather
using JLD2
using Unitful
drop = Base.Iterators.drop

push!(LOAD_PATH,joinpath(@__DIR__,"packages"))
using AuditoryModel
using AuditoryCoherence

# match these parameters better with the manual output I've been looking at
params = Dict(
  :condition => [:scales],
  :c_σ => 10.^linspace(-1,-0.1,5),  :τ_σ => linspace(0.0s,2.0s,3)[2:end],
  :c_m => 10.^linspace(-1,1,5),     :τ_m => linspace(50.0ms,500.0ms,3),
  :c_a => 10.^linspace(-1.0,1.0,5), :τ_a => linspace(100.0ms,10.0s,3),

  :W_m_σ => [2.0]
  :τ_x => linspace(100ms,300ms,2), :c_x => linspace(0.75,5,2),
  :c_n => [15], :τ_n => [1s],
)

function byparams(params)
  if length(params) == 1
    key,vals = first(params)
    DataFrame(;key => vals)
  else
    key,vals = first(params)
    others = byparams(drop(params,1))

    result = vcat((others for i in 1:length(vals))...)
    result[key] = repeat(vals,inner=nrow(others))
    result
  end
end
df = byparams(params)
categorical!(df,:condition)

in_ms(x) = Float64.(ustrip.(uconvert.(ms,x)))

open(joinpath(@__DIR__,"count_lengths_N.txt"),"w") do f
  println(f,"$(nrow(df))")
end

filename = joinpath(@__DIR__,"params.feather")
Feather.write(filename,df,transforms = Dict(
    "τ_σ" => in_ms,
    "τ_m" => in_ms,
    "τ_a" => in_ms,
    "τ_x" => in_ms,
    "τ_n" => in_ms
    "condition" => CategoricalArray
  ))

