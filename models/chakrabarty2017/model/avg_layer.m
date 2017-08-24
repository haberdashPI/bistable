function y=avg_layer(x)

for i=1:size(x,1)
    if(max(x(i,:))>0)
    y(i,:)=x(i,:)./max(x(i,:));
    else
     y(i,:)=-x(i,:)./max(x(i,:));  
    end
end