using Unitful
using Sounds
export alter_ab, sync_ab, buildup_aba, aba, ab, ideal_ab

set_default_samplerate!(8kHz)

function alter_ab(tone_len,ab_repeats,freq,delta)
  space = silence(tone_len,sample_rate=fs)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  space = silence(tone_len)

  ab_seq = silence(0s)
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

  ab_seq = silence(0s)
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

  aba_seq = silence(0s)
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

  aba_seq = silence(0s)
  for i = 1:repeats
    aba_seq = [aba_seq; a; gap; b; gap; a; gap; space; gap]
  end

  aba_seq
  # [a; space; a; space; aba_seq]
end


function ab(tone_len,spacing_len,offset_ratio,repeats,freq,delta,options...)
  @assert 0 <= offset_ratio <= 1
  space = silence(tone_len)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = :without_a in options ? silence(tone_len) : tone(a_freq,tone_len)
  b = :without_b in options ? silence(tone_len) : tone(b_freq,tone_len)
  a = :ramp in options ? ramp(a,25ms) : a
  b = :ramp in options ? ramp(b,25ms) : b
  ab = mix(a,[silence(offset_ratio*(tone_len+spacing_len));b])
  ab_len = 2(tone_len+spacing_len)
  ab = [ab; silence(ab_len - duration(ab))]

  ab_seq = silence(0s)
  for i = 1:repeats
    ab_seq = [ab_seq; ab]
  end

  ab_seq
end


function ideal_ab(spect,tone_len,spacing_len,offset_ratio,repeats,freq,
                  delta,options...)
  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a_freqi = indmin(abs.(freqs(spect) .- a_freq))
  b_freqi = indmin(abs.(freqs(spect) .- b_freq))

  a_times = 2(tone_len+spacing_len) .* (0:10)
  b_times = 2(tone_len+spacing_len) .* (0:10) .+
    offset_ratio*(tone_len+spacing_len)
  sp = zeros(ceil(Int,(maximum(b_times) + (tone_len+spacing_len)) / Δt(spect)),
             length(freqs(spect)))

  a_timei = find(mapslices(any,a_times .< times(spect,sp)' .<
                           (a_times .+ (tone_len)),[1]))
  b_timei = find(mapslices(any,b_times .< times(spect,sp)' .<
                           (b_times .+ (tone_len)),[1]))

  if !(:without_a in options) sp[a_timei,a_freqi] = 1 end
  if !(:without_b in options) sp[b_timei,b_freqi] = 1 end

  sp
end
