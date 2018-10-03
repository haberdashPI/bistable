using FileIO
using DataFrames
using Feather
using CategoricalArrays
using Unitful
using Unitful: ms, s, Hz, kHz
using Dates

drop = Base.Iterators.drop

in_ms(x) = Float64.(ustrip.(uconvert.(ms,x)))
in_Hz(x) = Float64.(ustrip.(uconvert.(Hz,x)))

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

unit_table = Dict(r"τ" => in_ms,r"Δt" => in_ms, r"^f$" => in_Hz)

function write_params(label,df)
  dir = joinpath(@__DIR__,"..","data","count_lengths","run_$(Date(now()))")
  isdir(dir) || mkdir(dir)
  open(joinpath(dir,"$(label)_N.txt"),"w") do f
    println(f,"$(size(df,1))")
  end

  filename = joinpath(dir,"$(label)_params.feather")

  for col in names(df)
    for (pattern,fn) in pairs(unit_table)
      if occursin(pattern,string(col))
        df[col] = fn.(df[col])
      end
    end
  end

  Feather.write(filename,df)
end

a_vals = [0.0; 10 .^ range(0.7,stop=1.75,length=4)]
m_vals = [0.0; 10 .^ range(1.25,stop=2,length=4)]

# this parameter set surveys all values for the 3st case. A second search will
# cover all three stimuli across any parameter sets that generate bistability
# for 3st

write_params("individual_levels",vcat(
  byparams(Dict(
    :Δt         => [240ms],
    :Δf         => [3,6,12],
    :f          => [500Hz],

    :f_c_x      => [3.0],    :f_τ_x     => [500ms],
    :f_c_σ      => [0.2],    :f_τ_σ     => [500ms],
    :f_c_a      => a_vals,   :f_τ_a     => [3s],
    :f_c_m      => m_vals,   :f_τ_m     => [350ms],
    :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

    :s_c_x      => [3.0],    :s_τ_x     => [500ms],
    :s_c_σ      => [0.0],    :s_τ_σ     => [500ms],
    :s_c_a      => [0.0],   :s_τ_a     => [3s],
    :s_c_m      => [0.0],   :s_τ_m     => [350ms],
    :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

    :t_c_x      => [3.0],    :t_τ_x     => [500ms],
    :t_c_σ      => [0.0],    :t_τ_σ     => [500ms],
    :t_c_a      => [0.0],  :t_τ_a     => [3s],
    :t_c_m      => [0.0],  :t_τ_m     => [350ms],
    :t_W_m_σ_t  => [7.0],  :t_W_m_σ_ϕ => [7.0],
    :t_W_m_σ_N  => [3.0],  :t_W_m_c    => [6.0]
  )),
  byparams(Dict(
    :Δt         => [240ms],
    :Δf         => [3,6,12],
    :f          => [500Hz],

    :f_c_x      => [3.0],    :f_τ_x     => [500ms],
    :f_c_σ      => [0.0],    :f_τ_σ     => [500ms],
    :f_c_a      => [0.0],   :f_τ_a     => [3s],
    :f_c_m      => [0.0],   :f_τ_m     => [350ms],
    :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

    :s_c_x      => [3.0],    :s_τ_x     => [500ms],
    :s_c_σ      => [0.2],    :s_τ_σ     => [500ms],
    :s_c_a      => a_vals,   :s_τ_a     => [3s],
    :s_c_m      => m_vals,   :s_τ_m     => [350ms],
    :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

    :t_c_x      => [3.0],    :t_τ_x     => [500ms],
    :t_c_σ      => [0.0],    :t_τ_σ     => [500ms],
    :t_c_a      => [0.0],   :t_τ_a     => [3s],
    :t_c_m      => [0.0],   :t_τ_m     => [350ms],
    :t_W_m_σ_t  => [7.0],  :t_W_m_σ_ϕ => [7.0],
    :t_W_m_σ_N  => [3.0],  :t_W_m_c    => [6.0]
  )),
  byparams(Dict(
    :Δt         => [240ms],
    :Δf         => [3,6,12],
    :f          => [500Hz],

    :f_c_x      => [3.0],    :f_τ_x     => [500ms],
    :f_c_σ      => [0.0],    :f_τ_σ     => [500ms],
    :f_c_a      => [0.0],   :f_τ_a     => [3s],
    :f_c_m      => [0.0],   :f_τ_m     => [350ms],
    :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

    :s_c_x      => [3.0],    :s_τ_x     => [500ms],
    :s_c_σ      => [0.0],    :s_τ_σ     => [500ms],
    :s_c_a      => [0.0],   :s_τ_a     => [3s],
    :s_c_m      => [0.0],   :s_τ_m     => [350ms],
    :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

    :t_c_x      => [3.0],    :t_τ_x     => [500ms],
    :t_c_σ      => [0.2],    :t_τ_σ     => [500ms],
    :t_c_a      => a_vals,   :t_τ_a     => [3s],
    :t_c_m      => m_vals,   :t_τ_m     => [350ms],
    :t_W_m_σ_t  => [7.0],  :t_W_m_σ_ϕ => [7.0],
    :t_W_m_σ_N  => [3.0],  :t_W_m_c    => [6.0]
  )),
 ))

# write_params("individual",vcat(
#   byparams(Dict(
#     :Δt         => [240ms],
#     :Δf         => [3],
#     :f          => [500Hz],

#     :f_c_x      => [3.0],    :f_τ_x     => [500ms],
#     :f_c_σ      => [0.2],    :f_τ_σ     => [500ms],
#     :f_c_a      => a_vals,   :f_τ_a     => [3s],
#     :f_c_m      => m_vals,   :f_τ_m     => [350ms],
#     :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

#     :s_c_x      => [3.0],    :s_τ_x     => [500ms],
#     :s_c_σ      => [0.0],    :s_τ_σ     => [500ms],
#     :s_c_a      => [0.0],   :s_τ_a     => [3s],
#     :s_c_m      => [0.0],   :s_τ_m     => [350ms],
#     :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

#     :t_c_x      => [3.0],    :t_τ_x     => [500ms],
#     :t_c_σ      => [0.0],    :t_τ_σ     => [500ms],
#     :t_c_a      => [0.0],   :t_τ_a     => [3s],
#     :t_c_m      => [0.0],   :t_τ_m     => [350ms],
#     :t_W_m_σ_t  => [8.0],    :t_W_m_σ_ϕ => [6.0],
#     :t_W_m_c    => [6.0]
#   )),
#   byparams(Dict(
#     :Δt         => [240ms],
#     :Δf         => [3],
#     :f          => [500Hz],

#     :f_c_x      => [3.0],    :f_τ_x     => [500ms],
#     :f_c_σ      => [0.0],    :f_τ_σ     => [500ms],
#     :f_c_a      => [0.0],   :f_τ_a     => [3s],
#     :f_c_m      => [0.0],   :f_τ_m     => [350ms],
#     :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

#     :s_c_x      => [3.0],    :s_τ_x     => [500ms],
#     :s_c_σ      => [0.2],    :s_τ_σ     => [500ms],
#     :s_c_a      => a_vals,   :s_τ_a     => [3s],
#     :s_c_m      => m_vals,   :s_τ_m     => [350ms],
#     :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

#     :t_c_x      => [3.0],    :t_τ_x     => [500ms],
#     :t_c_σ      => [0.0],    :t_τ_σ     => [500ms],
#     :t_c_a      => [0.0],   :t_τ_a     => [3s],
#     :t_c_m      => [0.0],   :t_τ_m     => [350ms],
#     :t_W_m_σ_t  => [8.0],    :t_W_m_σ_ϕ => [6.0],
#     :t_W_m_c    => [6.0]
#   )),
#   byparams(Dict(
#     :Δt         => [240ms],
#     :Δf         => [3],
#     :f          => [500Hz],

#     :f_c_x      => [3.0],    :f_τ_x     => [500ms],
#     :f_c_σ      => [0.0],    :f_τ_σ     => [500ms],
#     :f_c_a      => [0.0],   :f_τ_a     => [3s],
#     :f_c_m      => [0.0],   :f_τ_m     => [350ms],
#     :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

#     :s_c_x      => [3.0],    :s_τ_x     => [500ms],
#     :s_c_σ      => [0.0],    :s_τ_σ     => [500ms],
#     :s_c_a      => [0.0],   :s_τ_a     => [3s],
#     :s_c_m      => [0.0],   :s_τ_m     => [350ms],
#     :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

#     :t_c_x      => [3.0],    :t_τ_x     => [500ms],
#     :t_c_σ      => [0.2],    :t_τ_σ     => [500ms],
#     :t_c_a      => a_vals,   :t_τ_a     => [3s],
#     :t_c_m      => m_vals,   :t_τ_m     => [350ms],
#     :t_W_m_σ_t  => [8.0],    :t_W_m_σ_ϕ => [6.0],
#     :t_W_m_c    => [6.0]
#   )),
#  ))

# write_params("survey_interact",
# byparams(Dict(
#   :Δt         => [240ms],
#   :Δf         => [3],
#   :f          => [500Hz],

#   :f_c_x      => [3.0],    :f_τ_x     => [500ms],
#   :f_c_σ      => [0.2],    :f_τ_σ     => [500ms],
#   :f_c_a      => a_vals,   :f_τ_a     => [3s],
#   :f_c_m      => m_vals,   :f_τ_m     => [350ms],
#   :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

#   :s_c_x      => [3.0],    :s_τ_x     => [500ms],
#   :s_c_σ      => [0.2],    :s_τ_σ     => [500ms],
#   :s_c_a      => a_vals,   :s_τ_a     => [3s],
#   :s_c_m      => m_vals,   :s_τ_m     => [350ms],
#   :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

#   :t_c_x      => [3.0],    :t_τ_x     => [500ms],
#   :t_c_σ      => [0.2],    :t_τ_σ     => [500ms],
#   :t_c_a      => a_vals,   :t_τ_a     => [3s],
#   :t_c_m      => m_vals,   :t_τ_m     => [350ms],
#   :t_W_m_σ_t  => [8.0],    :t_W_m_σ_ϕ => [6.0],
#   :t_W_m_c    => [6.0]
#  ))
# )

