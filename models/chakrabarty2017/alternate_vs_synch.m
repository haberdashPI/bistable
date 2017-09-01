% clc;
% clear all
addpath('stimulus');
addpath('model');

loadload;

[~,hostname] = system('hostname')
if startsWith(hostname,'Claude')
  base_dir = '/Volumes/Miguel/Research/deb';
else
  base_dir = '/Volumes/Data/Little_Bistable_2017_08_15/deb';
end

aud_model = init_model(base_dir,1);

fs=8000;
tt = [1:0.06*fs]/fs;

deltas = [1 3 6 9 15];
freqs=[500 750 1000 1250 1500 1750 2000];

ab_repeats = 17;
tone_len = 0.06;

taus = [1,3];

a_dist = [];
b_dist = [];

for freq_i=1:length(freqs)

  % delta_i is for different Δf's between A and B
  for delta_i=1:length(deltas)
    for k = 1:2
      if k == 1 % alternating
        [base, a_reference, b_reference] = ...
            alternating_abab(tone_len,ab_repeats,fs,...
                             freqs(freq_i),deltas(delta_i));
      else % synchronous
        [base, a_reference, b_reference] = ...
            synchronous_abab(tone_len,ab_repeats,fs,...
                             freqs(freq_i),deltas(delta_i));
      end

      hebb_base = run_model(aud_model,taus,base);
      hebb_a_reference = run_model(aud_model,taus,a_reference);
      hebb_b_reference = run_model(aud_model,taus,b_reference);

      for tau = taus
        a_dist(delta_i,tau,freq_i,k) = ...
            hebb_dist(hebb_base{tau}(end-5:end,:),...
                      hebb_a_reference{tau}(end-5:end,:));
        b_dist(delta_i,tau,freq_i,k) = ...
            hebb_dist(hebb_base{tau}(end-5:end,:),...
                      hebb_b_reference{tau}(end-5:end,:));
      end
    end

    delta_i
  end

  freq_i
end

alter = mean(ab_resp(a_dist(:,:,:,1),b_dist(:,:,:,1),taus),1);
sync = mean(ab_resp(a_dist(:,:,:,2),b_dist(:,:,:,2),taus),1);

plot([alter;sync]')
legend('alternating','synchronous')

csvwrite('../../data/alter_v_sync.csv',[alter;sync;deltas]')
