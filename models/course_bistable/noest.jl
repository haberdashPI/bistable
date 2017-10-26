using Unitful: s,ms,Hz,kHz,ustrip

function simulate(input,W,dt)
  N = size(W,2)
  u = (1:N) + 0N
  S = reshape((1:N^2) + 1N,N,N)

  β = -2*(I - ones(N,N)) / (N-1)
  γ = 1e-3
  τ_u = 1ms
  τ_s = 1s
  τ_a = 1s
  φ = 1

  x = 1e-6rand(size(input,1),n_states)
  f(x) = max(zero(x),sqrt(x))
  for t in 2:size(x,1)
    x[t,:] = x[t-1,:]

    x[t,u] .+= (f.(input[t,:]*W .+ x[t-1,S]*β*x[t-1,u] .- γ*x[t-1,a]) .+
                -x[t,u]) ./ τ_u .* dt
    x[t,a] .+= (-x[t,a] + x[t,u]) ./ τ_a .* dt
    x[t,S] .+= (1 .- x[t,S] .- φ.*x[t,S].*x[t-1,u]') ./ τ_s .* dt
  end

  x
end

x = zeros(2,)
