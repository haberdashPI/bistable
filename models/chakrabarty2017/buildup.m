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

deltas = [1 3 6 9];
freqs=[500]; % 750 1000 1250 1500 1750 2000];
seq_dur=[500 1000 1500 2000 2500 3000 3500 4000 4500 ...
         5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000];

tone_len = 0.125;

taus = 1:3;

dist = [];

for freq_i=1:length(freqs)
  for delta_i=1:length(deltas)
    for seq_dur_i=1:length(seq_dur)
      % alternating stimulus
      [a_reference, b_reference] = ...
          buildup_aba(tone_len,ceil(seq_dur(seq_dur_i)/500),fs,...
                      freqs(freq_i),deltas(delta_i));

      hebb_a_reference = run_model(aud_model,taus,a_reference);
      hebb_b_reference = run_model(aud_model,taus,b_reference);

      cur_dist = 0
      for tau = taus
        cur_dist = cur_dist + ...
            hebb_dist(hebb_a_reference{tau}(end-5:end,:),...
                      hebb_b_reference{tau}(end-5:end,:));
      end

      dist(delta_i,seq_dur_i,freq_i) = cur_dist

      delta_i
    end

    freq_i
  end
end

mdist = mean(dist,3)
mdist = mdist(:,1) - mdist

t = 4;
sigma = 1;
resps = bsxfun(@gt,mdist,t+sigma*randn([size(mdist) 10000]));
prop_streaming = mean(resps,3);
prop_streaming = [prop_streaming; seq_dur];
csvwrite('../../data/buildup.csv',prop_streaming');
