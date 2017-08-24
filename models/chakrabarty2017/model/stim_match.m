function match = stim_match(dist,taus)
mu = mean(dist,3);
match = sum(mu(:,taus),2);
end
