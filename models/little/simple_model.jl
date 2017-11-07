using Unitful: s, ms, ustrip
# question: should we somehow limit the unit response???
# e.g. with:
sig(x) = 1/(1+exp(-10(x-0.5)))


function adaptmi_old(fn,x,τ_a,c_a,c_mi,shape=identity,Δt=1ms)
  y = similar(x)
  a = similar(x)

  y_t = zeros(size(x,2))
  a_t = zeros(size(x,2))

  c_mi = c_mi / (length(y_t)-1)
  dt_a = Δt / τ_a
  for t in 1:size(x,1)
    y_t = shape.(fn(x[t,:]) .- c_a.*a_t .- c_mi.*(sum(y_t) .- y_t))
    a_t .+= dt_a .* (y_t .- a_t)

    y[t,:] .= y_t
    a[t,:] = a_t
  end

  y,a
end

function adaptmi(fn,x::AbstractArray{T};
                 τ_a=1.5s,c_a=5,τ_mi=50ms,c_mi=10,
                 shape::Function=identity,Δt=1ms,
                 W_mi=ones(T,size(x,2),size(x,2)) - I) where T
  y = similar(x)
  a = similar(x)
  mi = similar(x)

  y_t = zeros(T,size(x,2))
  a_t = zeros(T,size(x,2))
  mi_t = zeros(T,size(x,2))

  c_mi = c_mi / (length(y_t)-1)
  dt_a = Δt / τ_a
  dt_mi = Δt / τ_mi
  for t in 1:size(x,1)
    y_t = shape.(fn(x[t,:]) .- c_a.*a_t .- c_mi.*mi_t)
    a_t .+= dt_a .* (y_t .- a_t)
    mi_t .+= dt_mi .* (W_mi*y_t .- mi_t)

    y[t,:] .= y_t
    a[t,:] = a_t
    mi[t,:] = mi_t
  end

  y,a,mi
end

function noise!(x,τ_ε,c_ε,Δt=1ms)
  ε_t = zeros(size(x,2))
  n = size(x,2)
  for t in 1:size(x,1)
    ε_t .+= -ε_t.*(Δt/τ_ε) + randn(n).*(c_ε*sqrt(2Δt/τ_ε))
    x[t,:] += ε_t
  end
  x
end
noise(x,τ_ε,c_ε) = noise!(copy(x),τ_ε,c_ε)
