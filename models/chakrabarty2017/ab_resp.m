function responses = ab_resp(a_dist,b_dist,taus)
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
