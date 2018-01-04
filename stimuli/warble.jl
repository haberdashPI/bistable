include("../models/little/audio_spect.jl")
setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=10)

function warble(freq,width,depth,scale,rate,len,itr=25,
                unit_len=min(3s,len),fade_len=1s)
  @assert fade_len < unit_len

  N = ceil(Int,unit_len / Δt(spect))
  y = 0.1ones(N,128)
  lower = freq * 2^(-width/2)
  upper = freq * 2^(width/2)
  ixs = find(lower .< freqs(spect,y) .< upper)
  r = linspace(0,2,length(ixs))

  y[:,ixs] .= (1 .- (r .- 1).^2)'

  t = ustrip.(indices(y,1)*Δt(spect))
  f = indices(y,2)*Δf(spect)

  shift(x) = x*0.5depth + (1-0.5depth)
  ym = y .* shift.(sin.(2π.*(scale.*f' .- rate.*t)))
  x = inv(spect,ym,iterations=itr)
  xs = @> sound(x) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2))

  if unit_len < len
    extender(x,n) = n <= 1 ? x : fadeto(x,extender(x,n-1),fade_len)
    extender(xs,ceil(Int,len / (unit_len-fade_len)))[0s .. len]
  else
    xs
  end
end

function noise_warble(freq,width,depth,scale,rate,len,itr=25,
                unit_len=min(3s,len),fade_len=1s)
  @assert fade_len < unit_len

  x = @> noise(unit_len) bandpass(freq * 2.0^(-width/2),freq * 2.0^(width/2))
  y = spect(x)

  t = ustrip.(indices(y,1)*Δt(spect))
  f = indices(y,2)*Δf(spect)

  shift(x) = x*0.5depth + (1-0.5depth)
  ym = y .* shift.(sin.(2π.*(scale.*f' .- rate.*t)))
  x = inv(spect,ym,iterations=itr)

  if unit_len < len
    extender(x,n) = n <= 1 ? x : fadeto(x,extender(x,n-1),fade_len)
    extender(sound(x),ceil(Int,len / (unit_len-fade_len)))[0s .. len]
  else
    x
  end
end
