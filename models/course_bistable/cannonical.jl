using Plots; plotlyjs()
using Unitful: s,ms,Hz,kHz,ustrip

baseline = 1e-4
on = 0.375

function simulate(input,W,C,dt)
  N = size(W,1)
  u = (1:N) + 0N
  a = (1:N) + 1N
  S = reshape((1:N^2) + 2N,N,N)
  n_states = 2N + N^2

  β = 2C / (N-1)
  γ = 1e-3
  τ_u = 30ms
  τ_s = τ_u*100
  τ_a = τ_u*100
  φ = 1

  x = zeros(size(input,1),n_states)
  x[1,u] = rand(length(u))
  x[:,S] = 1
  f(x) = x > zero(x) ? sqrt(x) : zero(x)
  for t in 2:size(x,1)
    x[t,:] = x[t-1,:]

    xhat = W*input[t,:] .+ (β.*x[t-1,S])*x[t-1,u] .- γ*x[t-1,a]
    x[t,u] .+= (f.(xhat) .+ -x[t-1,u]) ./ τ_u .* dt
    x[t,a] .+= (-x[t-1,a] + x[t-1,u]) ./ τ_a .* dt
    x[t,S] .+= (1 .- x[t,S] .- φ.*x[t-1,S].*x[t-1,u]') ./ τ_s .* dt
  end

  x
end
# note the inputs need to be similar to the two units
# because that is what leads to the multi-stable behavior
# in this model
off = 1.0 - 1e-2
W = [1.0 off
     off 1.0]
C = -[0 1
      1 0]

dt = 1.0ms
ts = (0s):dt:(5s)
input = baseline*ones(length(ts),2)
input[ts .% 0.5s .< 150ms,1] .= on
input[(ts .+ 0.25s) .% 0.5s .< 150ms,2] .= on


r1 = simulate(input,W,C,dt)


W = [1.0    1.0off
     0.8    0.8off
     1.0off 1.0
     0.8off 0.8]

C = -[0 0 1 1
      0 0 1 1
      1 1 0 0
      1 1 0 0]

r2 = simulate(input,W,C,dt)


W = [1.0    1.0off
     0.8    0.8off
     1.0off 1.0
ch     0.8off 0.8
     0.0    0.0
     0.0    0.0]

C = -2[0 0 1 1 0 0
      0 0 1 1 0 0
      1 1 0 0 0 0
      1 1 0 0 0 0
      0 0 0 0 0 0
      0 0 0 0 0 0]

r3 = simulate(input,W,C,dt)


W = [1.0    1.0off
     -0.8    -0.8off
     1.0off 1.0
     -0.8off -0.8
     0.0    0.0
     0.0    0.0]

C = -5[0 0 1 1 0 0
      0 0 1 1 0 0
      1 1 0 0 0 0
      1 1 0 0 0 0
      0 0 0 0 0 0
      0 0 0 0 0 0]

r4 = simulate(input,W,C,dt)
