using IterTools
using DataFrames
using ProgressMeter

include("rankin2015.jl")
using Rankin2015

function product_df(;kwds...)
  total = prod(length(k[2]) for k in kwds)
  df = DataFrame(;((k[1],similar(k[2],total)) for k in kwds)...)

  for (r,tuple) in enumerate(product(map(k -> k[2],kwds)...))
    for (j,k) in enumerate(kwds)
      df[r,k[1]] = tuple[j]
    end
  end
  return df
end

#########################################
## Parameter γ
df = product_df(
  DF = logspace(log10(1),log10(20),7),
  PR = logspace(log10(1),log10(20),7),
  repeat = 1:10,
  gamma = logspace(-2,-0.8,10)
)

sim_length = 240
df[:,:x] = 0.0
@showprogress 1 "Simulating [noise level]..." for r in 1:size(df,1)
  sim = Simulation(df[r,:DF],df[r,:PR],sim_length,
                   dt=1/200,params=BistableParams(γ=df[r,:gamma]))
  df[r,:x] = mean(responses(sim))
end

filename = "noise_level_"*Dates.format(now(),"yyyy-mm-dd_HH.MM")*".csv"
println("Wrote results to "*filename)
writetable(filename,df)


#########################################
## Parameter σ_p

df = product_df(
  DF = logspace(log10(1),log10(20),7),
  PR = logspace(log10(1),log10(20),7),
  repeat = 1:10,
  sigma_p = linspace(1,20,10)
)

sim_length = 240
df[:,:x] = 0.0
@showprogress 1 "Simulating [input spread]..." for r in 1:size(df,1)
  sim = Simulation(df[r,:DF],df[r,:PR],sim_length,
                   dt=1/200,params=BistableParams(σ_p=df[r,:sigma_p]))
  df[r,:x] = mean(responses(sim))
end
filename = "input_spread_"*Dates.format(now(),"yyyy-mm-dd_HH.MM")*".csv"
println("Wrote results to "*filename)
writetable(filename,df)

########################################
# Phase lengths
df = DataFrame()
n_repeats = 50
sim_length = 240
dt = 1/200
@showprogress 1 "Simulating [phase lengths]..." for r in 1:n_repeats
  sim = Simulation(5,8,sim_length,dt=dt)
  lengths = rlengths(responses(sim))
  df_0 = DataFrame(repeat = fill(r,length(lengths)),
                   lengths = lengths * dt,
                   index = indices(lengths,1))
  df = vcat(df,df_0)
end

filename = "phase_lengths_"*Dates.format(now(),"yyyy-mm-dd_HH.MM")*".csv"
println("Wrote results to "*filename)
writetable(filename,df)
