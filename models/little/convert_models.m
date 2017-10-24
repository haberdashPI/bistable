% converts models to julia readable format
base_dir = '/Volumes/Data/Little_Bistable_2017_08_15/deb';
addpath('model')
% aud_model = init_model(base_dir,1);

result_dir = '/Users/davidlittle/Data/';
file = [result_dir 'model.h5']
layer_1 = aud_model.layer_1.model;

h5create(file,'/layer1/W',size(layer_1.W{1}));
h5write(file,'/layer1/W', layer_1.W{1});

h5create(file,'/layer1/b',size(layer_1.biases{2}'));
h5write(file,'/layer1/b',layer_1.biases{2}');

h5create(file,'/layer1/bv',size(layer_1.biases{1}'));
h5write(file,'/layer1/bv',layer_1.biases{1}');

for tau_i = 1:8
  layer_2 = aud_model.layer_2{tau_i};

  h5create(file,['/layer2/tau' num2str(tau_i) '/W'],size(layer_2.vishid));
  h5write(file,['/layer2/tau' num2str(tau_i) '/W'],layer_2.vishid);

  h5create(file,['/layer2/tau' num2str(tau_i) '/Wpast'],size(layer_2.pasthid));
  h5write(file,['/layer2/tau' num2str(tau_i) '/Wpast'],layer_2.pasthid);

  h5create(file,['/layer2/tau' num2str(tau_i) '/b'],size(layer_2.hidbiases'));
  h5write(file,['/layer2/tau' num2str(tau_i) '/b'],layer_2.hidbiases');

  h5create(file,['/layer2/tau' num2str(tau_i) '/bv'],size(layer_2.visbiases'));
  h5write(file,['/layer2/tau' num2str(tau_i) '/bv'],layer_2.visbiases');
end
