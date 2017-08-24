function y = hebb_dist(a,b)
diff = a-b;
y = sqrt(trace(diff*diff'));