using FileIO
using JLD2
using DataFrames
using Feather

include("count_lengths.jl")
params = load("params.jld2","df")

# start by just looking the mean and max length
# of the two percepts across the various
# parameters

dir = joinpath("..","..","data","count_lengths","data")
all_rows = []
for file in readdir(dir)
  if ismatch(r"clbin$",file)
    rows = loadrows(joinpath(dir,file))
    if !any(x -> x.method ∉ 1:3,rows)
      push!(all_rows,DataFrame(pindex = map(x -> x.param_index,rows),
                               length = map(x -> x.length,rows),
                               stimulus = map(x -> x.stimulus,rows),
                               method = map(x -> method_labels[x.method],rows)))
    end
  end
end
df = vcat(all_rows...)

Feather.write(joinpath("..","..","data","count_lengths","all_rows.feather"),df)

r_params = DataFrame(
  c_sigma = params[:c_σ],
  tau_sigma = ustrip.(uconvert.(ms,params[:τ_σ])),
  c_m = params[:c_m],
  tau_m = ustrip.(uconvert.(ms,params[:τ_m])),
  c_a = params[:c_a],
  tau_a = ustrip.(uconvert.(ms,params[:τ_a])),
  W_m_sig = params[:W_m_σ],
  tau_x = ustrip.(uconvert.(ms,params[:τ_x])),
  c_x = params[:c_x],
  c_n = params[:c_n],
  τ_n = Float64.(ustrip.(uconvert.(ms,params[:τ_n]))),
)
Feather.write(joinpath("..","..","data","count_lengths","params.feather"),r_params)


