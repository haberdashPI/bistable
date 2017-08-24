function [base, a_reference, b_reference] = ...
    synchronous_abab(tone_len,ab_repeats,fs,freq,delta)

a_freq=freq;
b_freq=freq.*(2.^(delta/12));

a = tone(a_freq,tone_len,fs);
b = tone(b_freq,tone_len,fs);
space = silence(tone_len,fs);

ab_seq = [];
for i = 1:ab_repeats
  ab_seq = [ab_seq a+b space];
end

base = [ a space a space ab_seq ];
a_reference = [ base a space a space ];
b_reference = [ base b space b space ];
end
