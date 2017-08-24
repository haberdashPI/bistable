% freq_i is for different base F's
function responses = ab_stream_response(aud_model,freqs,deltas,stimfn,taus)
a_dist = [];
b_dist = [];

for freq_i=1:length(freqs)

  % delta_i is for different Δf's between A and B
  for delta_i=1:length(deltas)
    % create the stimulus
    [base, a_reference, b_reference] = stimfn(freqs(freq_i),deltas(delta_i));

    % run the model
    layer_1_base = run_model_layer(base,aud_model,1);
    layer_1_a_reference = run_model_layer(a_reference,aud_model,1);
    layer_1_b_reference = run_model_layer(b_reference,aud_model,1);

    for tau_i = taus
      layer_2_base = run_model_layer(layer_1_base,aud_model,2,tau_i);
      layer_2_a_reference = run_model_layer(layer_1_a_reference,aud_model,2,tau_i);
      layer_2_b_reference = run_model_layer(layer_1_b_reference,aud_model,2,tau_i);

      hebb_base = run_model_layer(layer_2_base,aud_model,3);
      hebb_a_reference = run_model_layer(layer_2_a_reference,aud_model,3);
      hebb_b_reference = run_model_layer(layer_2_b_reference,aud_model,3);

      a_dist(delta_i,tau_i,freq_i) = ...
          hebb_dist(hebb_base(end-5:end,:),hebb_a_reference(end-5:end,:));

      b_dist(delta_i,tau_i,freq_i) = ...
          hebb_dist(hebb_base(end-5:end,:),hebb_b_reference(end-5:end,:));

      tau_i
    end

    delta_i
  end

  freq_i
end

a_match = stim_match(a_dist,taus);
b_match = stim_match(b_dist,taus);

% correct=zeros(5,11);
% false=zeros(5,11);

t = 0.8;
sigma = 0.08;
N = 100;

a_resp = sum(bsxfun(@lt,a_match',t+sigma*randn(11,5,N)),3);
b_resp = sum(bsxfun(@lt,b_match',t+sigma*randn(11,5,N)),3);

responses = bsxfun(@dprime_simple,a_resp/(N+1),b_resp/(N+1));
end
