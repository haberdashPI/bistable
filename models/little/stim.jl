using Unitful
import Unitful: s, Hz

insamples(len::Quantity,fs) = floor(Int,ustrip(uconvert(s,len))*fs)
silence(tone_len,fs) = zeros(insamples(tone_len,fs))
function tone(freq::Quantity,len::Quantity,fs)
  t = (1:insamples(len,fs)) / fs
  sin.(2π * ustrip(uconvert(Hz,freq)) * t)
end

function alter_ab(tone_len,ab_repeats,fs,freq,delta)
  space = silence(tone_len,fs)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len,fs)
  b = tone(b_freq,tone_len,fs)
  space = silence(tone_len,fs)

  ab_seq = Float64[]
  for i = 1:ab_repeats
    ab_seq = [ab_seq; a; b]
  end

  base = [ a; space; a; space; ab_seq ]
  a_reference = [ base; a; space; a; space ]
  b_reference = [ base; b; space; b; space ]

  base, a_reference, b_reference
end


function sync_ab(tone_len,ab_repeats,fs,freq,delta)
  space = silence(tone_len,fs)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len,fs)
  b = tone(b_freq,tone_len,fs)
  space = silence(tone_len,fs)

  ab_seq = Float64[]
  for i = 1:ab_repeats
    ab_seq = [ab_seq; a+b; space]
  end

  base = [ a; space; a; space; ab_seq ]
  a_reference = [ base; a; space; a; space ]
  b_reference = [ base; b; space; b; space ]

  base, a_reference, b_reference
end

function buildup_aba(tone_len,repeats,fs,freq,delta)
  space = silence(tone_len,fs)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len,fs)
  b = tone(b_freq,tone_len,fs)
  space = silence(tone_len,fs)

  aba_seq = Float64[]
  for i = 1:repeats
    aba_seq = [aba_seq; a; b; a; space]
  end

  a_ref = [a; space; a; space; aba_seq; a; space; a; space]
  b_ref = [a; space; a; space; aba_seq; b; space; b; space]

  a_ref, b_ref
end


function aba(tone_len,repeats,fs,freq,delta)
  space = silence(tone_len,fs)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len,fs)
  b = tone(b_freq,tone_len,fs)
  space = silence(tone_len,fs)

  aba_seq = Float64[]
  for i = 1:repeats
    aba_seq = [aba_seq; a; b; a; space]
  end

  aba_seq
  # [a; space; a; space; aba_seq]
end
