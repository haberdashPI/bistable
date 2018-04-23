include("count_lengths.jl")

dir = joinpath("..","..","data","count_lengths")
isdir(dir) || mkdir(dir)

N = 1

methods = Dict(
  :peaks => x -> source_count_by_peaks(x,window=1s,delta=0.25s,buildup=1s),
  :threshold => x -> source_count_by_threshold(x,window=1s,delta=0.25s,
                                               cuttoff=2cycoct,buildup=1s)
)

params = Dict(
    :c_m=>5,:τ_m=>350ms,:W_m_σ=>1.0,
    :c_d=>1.8,:τ_d=>500ms,
    :c_e=>1.5,:τ_e=>300ms,
    :c_a=>30,:τ_a=>3s,:shape_y => x -> max(0,x),
    :α=>3.0,

    :τ_σ => 500ms,:c_σ => 0.3,

    :scales=>cycoct.*round.(2.0.^linspace(-1,2,9),1)
)

count_lengths("test",dir,N,methods,params)
# x = ab(120ms,120ms,1,50,500Hz,6) |> normpower |> amplify(-10dB)
# rows = count_lengths_helper(x,methods,params)
