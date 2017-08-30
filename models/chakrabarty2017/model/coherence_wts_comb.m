function y = coherence_wts_comb(data)

thresh = 0.9;

C = zeros(size(data,2));
n = size(C,1);
y = zeros(size(data));
for j=1:size(data,1)-1
  k = data(j,:) .* sign(data(j,:)-thresh);
  C = C + k' * k;
  C(1:(n+1):n*n) = 0; % clear diagonal

  y(j,:) = 2*avg_layer((C*data(j,:)')');
end
