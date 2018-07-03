using FileIO
using DataFrames
using JLD2
using Feather
using CategoricalArrays
using Unitful
drop = Base.Iterators.drop

push!(LOAD_PATH,joinpath(@__DIR__,"..","packages"))
using AuditoryModel
using AuditoryCoherence

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

df = vcat(
byparams(Dict(
  :delta_t    => [240ms],                          :delta_f     => [2,6,12],
  :standard_f => [500Hz],                          :condition   => [:scales],
  :c_x        => [3.0],                            :τ_x         => [500ms],
  :c_σ        => linspace(0,1,6),                  :τ_σ         => [500ms],
  :c_a        => [0.0;10.^linspace(0.75,1.75,5)],  :τ_a         => [3s],
  :c_m        => [0.0;10.^linspace(1.25,2,5)],     :τ_m         => [350ms],
  :W_m_σ      => [15.0],                           :W_m_c       => [6.0]
 )),
byparams(Dict(
  :delta_t    => [240ms],                          :delta_f     => [2,6,12],
  :standard_f => [500Hz],                          :condition   => [:freqs],
  :c_x        => [3.0],                            :τ_x         => [500ms],
  :c_σ        => [0.2],                            :τ_σ         => [500ms],
  :c_a        => [0.0;10.^linspace(0.75,1.75,5)],  :τ_a         => [3s],
  :c_m        => [0.0;10.^linspace(1.25,2,5)],     :τ_m         => [350ms],
  :W_m_σ      => linspace(5.0,100.0,6),            :W_m_c       => [6.0]
 )))
# categorical!(df,:condition)

in_ms(x) = Float64.(ustrip.(uconvert.(ms,x)))
in_Hz(x) = Float64.(ustrip.(uconvert.(Hz,x)))

open(joinpath(@__DIR__,"count_lengths_N.txt"),"w") do f
  println(f,"$(nrow(df))")
end

filename = joinpath(@__DIR__,"params_$(Date(now())).jld2")
save(filename,"params",df)

filename = joinpath(@__DIR__,"params_$(Date(now())).feather")

Feather.write(filename,df,transforms = Dict{String,Function}(
  "τ_x" => x -> in_ms.(x),
  "τ_σ" => x -> in_ms.(x),
  "τ_a" => x -> in_ms.(x),
  "τ_m" => x -> in_ms.(x),
  "delta_t" => x -> in_ms.(x),
  "standard_f" => x ->in_Hz.(x),
  "condition" => x -> string.(x)
 ))


