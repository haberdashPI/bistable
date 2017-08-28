function resp = calc_response_gen_mod(x,model,locs)
ncases = size(x,1);
ndims = size(x,2);

nt = model.nt;
batch_indices=nt+1:size(x,1);
nc=length(batch_indices);
data = single(x(batch_indices,:));
past = zeros(nc,nt*ndims,'single');

if nargin < 3
  locs = 1:(nt*ndims);
end

for hh=nt:-1:1 %note reverse order
  past(:,ndims*(nt-hh)+1:ndims*(nt-hh+1)) = x(batch_indices-hh,:);
end

resp=zeros(nc,model.numhid);

% combine responses across the mixture components
for cc=1:model.numcomp
  % why calculate these? they're not used right now....
  % v_biases = repmat(model.visbiases(cc,locs),nc,1) + ...
  %     past*model.pastvis(locs,locs,cc);
  h_biases = repmat(model.hidbiases(cc,:),nc,1)+ ...
      past*model.pasthid(locs,:,cc);

  resp=resp+(data*model.vishid(locs,:,cc) + h_biases);
end
