using Unitful: s,ms,Hz,kHz,ustrip

function generate_input(dt,stim_len,tone_len)
  input = zeros(floor(Int,stim_len/dt),2)

  Λ = 1/6      # plateau fraction
  α_1 = 0.015  # impulse rise time
  α_2 = 0.0825 # plateau stay time
  response(t) = (exp(2)/α_1^2) * t^2*exp(-2t / α_1) +
    (Λ*exp(2)/α_2^2) * t^2*exp(-2t / α_2)

  # the a's
  at = 0s
  t = 0s:dt:tone_len
  while at+tone_len < stim_len
    from = floor(Int,at / dt) + 1
    to = floor(Int,(at + tone_len) / dt) + 1
    input[from:to,1] = response.(ustrip.(uconvert.(s,(1:(to-from+1))*dt)))

    at += 2tone_len
  end

  # the b's
  at = tone_len
  while at+tone_len < stim_len
    from = floor(Int,at / dt) + 1
    to = floor(Int,(at + tone_len) / dt) + 1
    input[from:to,2] = response.(ustrip.(uconvert.(s,(1:(to-from+1))*dt)))

    at += 4tone_len
  end

  input
end

function simulate(input,dt,W)
  N = size(W,2)

  a = (1:N) + 0N
  n = (1:N) + 1N
  e = (1:N) + 2N
  r = (1:N) + 3N
  n_states = 4N

  τ_a = 1.4s
  τ_r = 10.0ms
  τ_n = 100.0ms
  τ_e = 150.0ms
  c_a = 1
  c_mi = 0.6
  c_n = 0.01

  M = c_mi*(I - ones(N,N)) / (N-1)

  x = 1e-6rand(size(input,1),n_states)

  σ(x) = 1 / (1 + exp(12(-x + 0.2)))
  for i in 2:size(x,1)
    x[i,:] = x[i-1,:]
    for k in 1:3
      x[i,r[k]] += (-x[i-1,r[k]] +
                    σ(input[i,:]'*W[:,k] + # input
                      x[i-1,e[k]] + # excitation
                      x[i-1,r]'*M[:,k] + # mutual inhibition
                      -c_a*x[i-1,a[k]] + # adaptation
                      x[i-1,n[k]])/ # noise
                    τ_r * dt)
      x[i,e[k]] += (-x[i-1,e[k]] + x[i-1,r[k]])/τ_e * dt
      x[i,a[k]] += (-x[i-1,a[k]] + x[i-1,e[k]])/τ_a * dt
      x[i,n[k]] += -x[i-1,n[k]] / τ_n * dt +
        c_n * sqrt(2/τ_n) * sqrt(dt) * randn()
    end
  end

  x
end

W = [1.0 0.0
     0.0 1.0
     0.5 0.5]'
dt = 1.0ms
t = 0s:dt:10s
x = simulate(generate_input(dt,50.0s,125.0ms),dt,W)
plot(ustrip(t),x[1:length(t),(1:3) + 3*3])
