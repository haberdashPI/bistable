% TODO: run and compare the results of the two versions

% clc;
% clear all
addpath('stimulus');
addpath('model');

loadload;

base_dir = '/Volumes/Data/Little_Bistable_2017_08_15/deb';
aud_model = init_model(base_dir,1);

fs=8000;
tt = [1:0.06*fs]/fs;

delta_f = [1 3 6 9 15];
freq_a=[500 750 1000 1250 1500 1750 2000];

ab_repeats = 17;
tone_len = 0.06;

taus = [1,3];

stimfn = @(f,d)alternating_abab(tone_len,ab_repeats,fs,f,d);
a_responses = ab_stream_response(aud_model,freq_a,delta_f,stimfn,taus);

stimfn = @(f,d)synchronous_abab(tone_len,ab_repeats,fs,f,d);
s_responses = ab_stream_response(aud_model,freq_a,delta_f,stimfn,taus);

alter = mean(a_responses,1);
sync = mean(s_responses,1);

plot([alter;sync]')
legend('atlernating','synchronous')
