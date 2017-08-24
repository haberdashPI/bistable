function out = run_model(x,model,layer,layer2_i)

if layer == 1
  spect = wav2aud(x,model.spect_paras);
  out = hidden_act(spect,model.layer_1.model);
  out = inertia_delta(out,model.delta_intertia);
  out = out(:,model.loc);
elseif layer == 2
  if model.preloaded
    layer_2i = model.layer_2{layer2_i};
  else
    layer_2i = load([model.layer_2_dir 'mod_ep' int2str(layer2_i) '.mat']);
  end

  filtered = rate_filter(x,layer2_i,size(x,2));
  layer_2_locs = loc_find(model.loc,layer2_i);
  out = calc_response_gen_mod(filtered,layer_2i,layer_2_locs);
  out = avg_layer(out);

  if ~model.preloaded
    clear layer_2i;
  end
else % layer == 3
  out = coherence_wts_comb(x,0);
end
