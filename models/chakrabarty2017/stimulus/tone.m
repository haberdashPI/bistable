function y=tone(freq,len,fs)
% TODO: should the tones have ramps?
t = (1:floor(len*fs)) / fs;
y = sin(2*pi*freq*t);