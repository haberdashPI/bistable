% freq_i is for different base F's
function responses = run_model(aud_model,taus,x)
% run the model
layer_1 = run_model_layer(x,aud_model,1);
for tau_i = taus
  layer_2 = run_model_layer(layer_1,aud_model,2,tau_i);
  responses{tau_i} = run_model_layer(layer_2,aud_model,3);

  tau_i
end

end
