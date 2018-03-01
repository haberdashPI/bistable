using Sounds
using Missings

function aba(;freq=nothing,df=nothing,tone_length=nothing,
             spacing=nothing,repeat=nothing)
  a = @> tone(freq,tone_length) normpower ramp
  if ismissing(df)
    b = silence(tone_length)
  else
    b = @> tone(freq*2^(df/12),tone_length) normpower ramp
  end

  space = silence(spacing)
  unit = [a;space;b;space;a;space;silence(duration([a;space]))]
  vcat(Iterators.repeated(unit,repeat)...)
end

function aba_pulse(;freq=nothing,pulse=nothing,pulse_spacing=nothing,
             df=nothing,tone_length=nothing,
             spacing=nothing,repeat=nothing)
  function pulsetone(freq,tone_length)
    pulse_length = (tone_length - pulse_spacing*(pulse-1)) / pulse
    ptone = @> tone(freq,pulse_length) ramp(pulse_spacing/2)
    pbreak = silence(pulse_spacing)
    vcat(Iterators.repeated([ptone;pbreak],pulse-1)...,ptone)
  end

  a = @> pulsetone(freq,tone_length) normpower ramp
  if ismissing(df)
    b = silence(tone_length)
  else
    b = @> pulsetone(freq*2^(df/12),tone_length) normpower ramp
  end

  space = silence(spacing)
  unit = [a;space;b;space;a;space;silence(duration([a;space]))]
end

function aba_band(;freq=nothing,width=nothing,df=nothing,tone_length=nothing,
                 spacing=nothing,repeat=nothing)
  a = @> begin
    noise(tone_length)
    bandpass(freq*2.0^(-width/2),freq*2.0^(width/2))
    normpower
    ramp
  end
  if ismissing(df)
    b = silence(tone_length)
  else
    b = @> begin
      noise(tone_length)
      bandpass(@show(freq*2.0^(df/12-width/2)),@show(freq*2.0^(df/12+width/2)))
      normpower
      ramp
    end
  end

  space = silence(50ms)
  unit = [a;space;b;space;a;space;silence(duration([a;space]))]
  vcat(Iterators.repeated(unit,repeat)...)
end
