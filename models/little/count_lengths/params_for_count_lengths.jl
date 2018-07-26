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
  :delta_t    => [240ms],                         :delta_f   => [0.5,3,12],
  :standard_f => [500Hz],                         :condition => [:scales_track],

  :s_c_x      => [3.0],                           :s_τ_x     => [500ms],
  :s_c_σ      => [0.2],                           :s_τ_σ     => [500ms],
  :s_c_a      => [0.0;10.^linspace(0.75,1.75,5)], :s_τ_a     => [3s],
  :s_c_m      => [0.0;10.^linspace(1.25,2,5)],    :s_τ_m     => [350ms],
  :s_W_m_σ    => [15.0],                          :s_W_m_c   => [6.0],

  :t_c_x      => [3.0],                           :t_τ_x     => [500ms],
  :t_c_σ      => [0.2],                           :t_τ_σ     => [500ms],
  :t_c_a      => [0.0;10.^linspace(0.75,1.75,5)], :t_τ_a     => [3s],
  :t_c_m      => [0.0;10.^linspace(1.25,2,5)],    :t_τ_m     => [350ms],
  :t_W_m_σ    => [15.0],                          :t_W_m_c   => [6.0],

  :θ          => [1.3,1.6,1.75,2.1,2.4]
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

