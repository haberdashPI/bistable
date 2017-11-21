using Unitful
using TimedSound

setup_sound(sample_rate=8kHz)

function alter_ab(tone_len,ab_repeats,freq,delta)
  space = silence(tone_len,sample_rate=fs)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  space = silence(tone_len)

  ab_seq = Float64[]
  for i = 1:ab_repeats
    ab_seq = [ab_seq; a; b]
  end

  base = [ a; space; a; space; ab_seq ]
  a_reference = [ base; a; space; a; space ]
  b_reference = [ base; b; space; b; space ]

  base, a_reference, b_reference
end


function sync_ab(tone_len,ab_repeats,freq,delta)
  space = silence(tone_len)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  space = silence(tone_len)

  ab_seq = Float64[]
  for i = 1:ab_repeats
    ab_seq = [ab_seq; a+b; space]
  end

  base = [ a; space; a; space; ab_seq ]
  a_reference = [ base; a; space; a; space ]
  b_reference = [ base; b; space; b; space ]

  base, a_reference, b_reference
end

function buildup_aba(tone_len,repeats,freq,delta)
  space = silence(tone_len)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  space = silence(tone_len)

  aba_seq = Float64[]
  for i = 1:repeats
    aba_seq = [aba_seq; a; b; a; space]
  end

  a_ref = [a; space; a; space; aba_seq; a; space; a; space]
  b_ref = [a; space; a; space; aba_seq; b; space; b; space]

  a_ref, b_ref
end


function aba(tone_len,repeats,freq,delta)
  space = silence(tone_len)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  space = silence(tone_len)

  aba_seq = Float64[]
  for i = 1:repeats
    aba_seq = [aba_seq; a; b; a; space]
  end

  aba_seq
  # [a; space; a; space; aba_seq]
end

function aba(tone_len,gap_len,repeats,freq,delta)
  space = silence(tone_len)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  gap = silence(gap_len)
  space = silence(tone_len)

  aba_seq = Float64[]
  for i = 1:repeats
    aba_seq = [aba_seq; a; gap; b; gap; a; gap; space; gap]
  end

  aba_seq
  # [a; space; a; space; aba_seq]
end
