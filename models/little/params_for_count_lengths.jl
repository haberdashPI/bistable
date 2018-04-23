using FileIO
using Feather
import Base.Iterators: drop

push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence

params = Dict(
  :c_σ => 10.^linspace(-1,-0.1,2),
  :τ_σ => linspace(0.0s,2.0s,2),
  :c_m => 10.^linspace(-1,1,2),
  :τ_m => linspace(50.0ms,500.0ms,2),
  :c_a => 10.^linspace(-1.0,1.0,2),
  :τ_a => linspace(100.0ms,10.0s,2),

  :W_m_σ => 10.^linspace(-1.0,1.0,2),
  :α => linspace(0.75,5,2),
  :c_e => 10.^linspace(-0.1,1,2),
  :τ_e => linspace(50.0ms,500.0ms,2),
  :c_d_diff => linspace(0.1,3,2),
  :τ_d_diff => linspace(50.0ms,100ms,2)
)

function withentry(dict,key,val)
  result = copy(dict)
  result[key] = val
  result
end

function byparams(params)
  if isempty(params)
    [Dict{Symbol,Number}()]
  else
    key,vals = first(params)
    others = byparams(drop(params,1))

    (withentry(dict,key,val) for val in vals for dict in others)
  end
end

Feater.write("params.feather",byparams(params))
