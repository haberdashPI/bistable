using FileIO
using JLD2
import Base.Iterators: drop

push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence

params = Dict(
  :c_σ => 10.^linspace(-1,-0.1,10),  :τ_σ => linspace(0.0s,2.0s,3)[2:end],
  :c_m => 10.^linspace(-1,1,10),     :τ_m => linspace(50.0ms,500.0ms,3),
  :c_a => 10.^linspace(-1.0,1.0,10), :τ_a => linspace(100.0ms,10.0s,3),

  :W_m_σ => 10.^linspace(-1.0,2.0,5),
  :τ_x => linspace(100ms,300ms,3), :c_x => linspace(0.75,5,3),
  :c_n => [15], :τ_n => [1s],
  :scale_start => [-1], :scale_stop => [2], :scale_N => [12]
)

function withentry(dict,key,val)
  result = copy(dict)
  result[key] = val
  result
end

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
save("params.jld2",:df => df)
