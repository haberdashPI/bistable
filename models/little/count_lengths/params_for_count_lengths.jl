using FileIO
using DataFrames
using JLD2
using Unitful
drop = Base.Iterators.drop

push!(LOAD_PATH,joinpath(@__DIR__,"..","packages"))
using AuditoryModel
using AuditoryCoherence

params = Dict(
  :delta_t    => [240ms],                          :delta_f     => [2,6,12],
  :standard_f => [500Hz],                          :condition   => [:scales],
  :c_x        => [3.0],                            :τ_x         => [500ms],
  :c_σ        => linspace(0,1,6),                  :τ_σ         => [500ms],
  :c_a        => [0.0;10.^linspace(0.75,1.75,5)],  :τ_a         => [10s],
  :c_m        => [0.0;10.^linspace(1.25,2,5)],     :τ_m         => [350ms],
  :W_m_σ      => [15.0],                           :W_m_c       => [6.0]
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
# categorical!(df,:condition)

in_ms(x) = Float64.(ustrip.(uconvert.(ms,x)))
in_Hz(x) = Float64.(ustrip.(uconvert.(Hz,x)))

open(joinpath(@__DIR__,"count_lengths_N.txt"),"w") do f
  println(f,"$(nrow(df))")
end

filename = joinpath(@__DIR__,"params_$(Date(now())).jld2")

save(filename,"params",df)

