using Sounds

freq = 500Hz
width = 4/12
a = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
freq = 500Hz*2^(6/12)
b = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
space = silence(50ms)
aba = [a;space;b;space;a;space;silence(duration([a;space]))]
scene = vcat(repeated(aba,10)...)
save("wide1_center.wav",scene)

freq = 500Hz
width = 2/12
a = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
freq = 500Hz*2^(6/12)
b = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
space = silence(50ms)
aba = [a;space;b;space;a;space;silence(duration([a;space]))]
scene = vcat(repeated(aba,10)...)
save("wide2_center.wav",scene)

freq = 500Hz
width = 1/12
a = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
freq = 500Hz*2^(6/12)
b = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
space = silence(50ms)
aba = [a;space;b;space;a;space;silence(duration([a;space]))]
scene = vcat(repeated(aba,10)...)
save("wide3_center.wav",scene)

freq = 500Hz
width = 4/12
a = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
freq = 500Hz*2^(10/12)
b = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
space = silence(50ms)
aba = [a;space;b;space;a;space;silence(duration([a;space]))]
scene = vcat(repeated(aba,10)...)
save("wide1_edge.wav",scene)

freq = 500Hz
width = 2/12
a = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
freq = 500Hz*2^(7/12)
b = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
space = silence(50ms)
aba = [a;space;b;space;a;space;silence(duration([a;space]))]
scene = vcat(repeated(aba,10)...)
save("wide2_edge.wav",scene)

freq = 500Hz
width = 1/12
a = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
freq = 500Hz*2^(6.25/12)
b = noise(95ms) |> bandpass(freq*2.0^(-width/2),freq*2.0^(width/2)) |>
  normpower |> ramp
space = silence(50ms)
aba = [a;space;b;space;a;space;silence(duration([a;space]))]
scene = vcat(repeated(aba,10)...)
save("wide3_edge.wav",scene)
