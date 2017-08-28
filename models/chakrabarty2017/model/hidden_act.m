function y=hidden_act(x,R)
steps = 3;
% the -1 is for consistentcy with the old implementation,
% if I remove this I can replce the below lines with
% rate_filter.
n = floor(size(x,1)/steps)-1;
x = reshape(x(1:steps*n,:)',[],n)';

x = bsxfun(@minus,x,mean(x,2));
x = bsxfun(@rdivide,x,std(x,[],2));

y=bsxfun(@plus,x*R.W{1}, R.biases{2}');
