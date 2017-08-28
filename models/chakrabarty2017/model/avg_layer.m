function y=avg_layer(x)
y = x;
% normalize by maximum
row_max = max(x,[],2);
y = y ./ row_max;

% keep sign of final result consistent with input
positive_rows = row_max > 0;
y(~positive_rows,:) = -y(~positive_rows,:);
