function y=aba(a_freq,b_freq,tone_len,silence_len,repeat,fs)

a = tone(a_freq,tone_len,fs);
b = tone(b_freq,tone_len,fs);
space = silence(silence_len,fs);
unit = [a space b space a];

spacing = 4*floor(tone_len*fs);
unit_len = length(unit);
y = zeros(1,spacing * repeat);
for i = 1:repeat
  from = (i-1)*spacing + 1;
  to = (from + unit_len - 1);
  y(from:to) = unit;
end
