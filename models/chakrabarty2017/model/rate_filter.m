function y=rate_filter(x,hist)

if(hist==1)
  y = x;
else
  n=floor(size(x,1)/hist);
  y=reshape(x(1:hist*n,:)',[],n)';
end
