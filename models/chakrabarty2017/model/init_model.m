function model = init_model(base_dir,preload_layer2)

model = struct();
model.spect_paras = [10 8 -2 -1];
model.layer_1=load([base_dir '/model.mat']);
model.delta_intertia = 4;

if ~(nargin == 1) && preload_layer2
  model.preloaded = 1;
  model.layer_2 = cell(1,8);

  layer_2_dir = [base_dir '/generic_mod_newest/'];
  dir=[base_dir '/generic_mod_newest/'];
  for i = 1:8;
    disp(['Loading layer 2, tau ' num2str(i)]);
    model.layer_2{i} = load([layer_2_dir 'mod_ep' int2str(i) '.mat']);
  end
else
  model.preloaded = 0;
  model.layer_2_dir = [base_dir '/generic_mod_newest/'];
end

model.loc=1:350;
