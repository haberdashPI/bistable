include("count_lengths.jl")

dir = joinpath("..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

N = 1
M = 1

methods = Dict(
  # :peaks => x -> source_count_by_peaks(x,window=1s,delta=0.25s,buildup=1s),
  :threshold => x -> source_count_by_threshold(x,window=1s,delta=0.25s,
                                               cutoff=2cycoct,buildup=1s)
)

params = Dict(
    :τ_x => 300ms, :c_x=> 3.0,
    :τ_n => 1s,    :c_n => 15,

    :τ_m => 350ms, :c_m => 30,  :W_m_σ=>10.0,
    :τ_a => 3s,    :c_a => 10, 
    :τ_σ => 500ms, :c_σ => 0.3,

    :scale_start => -1, :scale_stop => 2, :scale_N => 12
)

count_lengths("test",dir,N,M,methods,params)
# x = ab(120ms,120ms,1,50,500Hz,6) |> normpower |> amplify(-10dB)
# rows = count_lengths_helper(x,methods,params)
