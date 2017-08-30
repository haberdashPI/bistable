function y=avg_layer(x)
y = x;
% normalize by maximum
row_max = max(x,[],2);
y = bsxfun(@rdivide,y,row_max);

% keep sign of final result consistent with input
y(row_max <= 0,:) = -y(row_max <= 0,:);
