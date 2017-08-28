function y = inertia_delta(x,hist)
y_hist = x(1:hist,:);
high = y_hist > 30;
y_hist(high) = y_hist(high)*5;
y_hist(~high) = y_hist(~high)/5;

y = x;
y(1:hist,:) = y_hist;
