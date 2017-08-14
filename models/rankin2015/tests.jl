include("rankin2015.jl")
using Rankin2015

using IterTools
using DataFrames
using ProgressMeter

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
                   length = lengths,
                   index = indices(lengths,1))
  df = vcat(df,df_0)
end


#########################################
# van noorden for intermittant
df = product_df(
  DF = logspace(log10(1),log10(20),7),
  PR = logspace(log10(1),log10(20),7),
  repeat = 1:10,
)

sim_length = 240
df[:,:x] = 0.0
@showprogress 1 "Simulating [van Noorden]..." for r in 1:size(df,1)
  inputs = create_int_inputs(3,df[r,:DF],df[r,:PR],sim_length,dt=dt)
  sim = Simulation(df[r,:DF],df[r,:PR],sim_length,dt=1/200)
  df[r,:x] = mean(int_responses(3,df[r,:PR],sim))
end

filename = "../../int_van_noorden_"*Dates.format(now(),"yyyy-mm-dd_HH.MM")*".csv"
println("Wrote results to "*filename)
writetable(filename,df)


########################################
# Intermittant stimulus, continuous phase lengths
df = DataFrame()
n_repeats = 50
sim_length = 240
dt = 1/200
@showprogress 1 "Simulating [phase lengths]..." for r in 1:n_repeats
  inputs = create_int_inputs(3,5,8,sim_length,dt=dt)
  sim = Simulation(5,8,sim_length,dt=dt,inputs=inputs)
  lengths = rlengths(responses(sim))
  df_0 = DataFrame(repeat = fill(r,length(lengths)),
                   length = lengths,
                   index = indices(lengths,1))
  df = vcat(df,df_0)
end

filename = ("../../data/int_cont_phase_lengths_"*
              Dates.format(now(),"yyyy-mm-dd_HH.MM")*".csv")
println("Wrote results to "*filename)
writetable(filename,df)


########################################
# Intermittant phase lengths
df = DataFrame()
n_repeats = 50
sim_length = 240
dt = 1/200
@showprogress 1 "Simulating [phase lengths]..." for r in 1:n_repeats
  inputs = create_int_inputs(3,5,8,sim_length,dt=dt)
  sim = Simulation(5,8,sim_length,dt=dt,inputs=inputs)
  lengths = rlengths(int_responses(3,8,sim))
  df_0 = DataFrame(repeat = fill(r,length(lengths)),
                   length = lengths,
                   index = indices(lengths,1))
  df = vcat(df,df_0)
end

filename = ("../../data/int_phase_lengths_"*
              Dates.format(now(),"yyyy-mm-dd_HH.MM")*".csv")
println("Wrote results to "*filename)
writetable(filename,df)

########################################
## Intermittent rate graph
df = DataFrame()
sim_length = 20
dt = 1/200
inputs = create_int_inputs(3,5,8,sim_length,dt=dt)
t,u,If,rs = collect(Simulation(5,8,sim_length,dt=dt,inputs=inputs))

marks = sum(If .> 0.5,2) .> 0

using Plots; gr()

plot(t,u[:,Rankin2015.r],label=["A" "AB" "B"])
scatter!(t[find(marks)],zeros(sum(marks)),
         mswidth=0,mcolor="grey",label="input")

filename = ("../../plots/int_responses_"*
              Dates.format(now(),"yyyy-mm-dd")*".pdf")
println("Wrote plot to "*filename)

savefig(filename)
