function  x=create_tone(tone_a,tone_b,dur)

x=[];

for i=1:dur
  %  x=[x sil1 tone_a sil2 sil1 tone_a sil2];
     x=[x tone_a tone_b];
end