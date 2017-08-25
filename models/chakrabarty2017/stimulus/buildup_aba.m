function [a_ref,b_ref] = buildup_aba(tone_len,repeats,fs,freq,delta)
space = silence(tone_len,fs);

a_freq=freq;
b_freq=freq.*(2.^(delta/12));

a = tone(a_freq,tone_len,fs);
b = tone(b_freq,tone_len,fs);
space = silence(tone_len,fs);

aba_seq = [];
for i = 1:repeats
  aba_seq = [aba_seq a b a space];
end

a_ref = [ a space a space aba_seq a space a space ];
b_ref = [ a space a space aba_seq b space b space ];
end