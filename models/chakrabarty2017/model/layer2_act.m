function resp = layer2_act(x,model,locs)
if nargin < 3
  locs = 1:(nt*size(x,2));
end

resp=zeros(size(x,1)-1,model.numhid);

% combine responses across the mixture components
for cc=1:model.numcomp
  % why calculate these? they're not used right now....
  % v_biases = repmat(model.visbiases(cc,locs),size(x,1)-1,1) + ...
  %     past*model.pastvis(locs,locs,cc);
  h_biases = bsxfun(@plus,model.hidbiases(cc,:),...
                    x(1:end-1,:)*model.pasthid(locs,:,cc));

  resp=resp+(x(2:end,:)*model.vishid(locs,:,cc) + h_biases);
end
