function [resp_data]= inertia_delta(data,hist)

resp_data=data;
for n=1:350
  temp=resp_data(:,n);
  temp1=temp(1:hist);

  for i=1:length(temp1)
    if(temp(i)>30)
      temp(i)=temp(i)*5;
    else
      temp(i)=temp(i)/5;
    end
  end

  resp_data(:,n)=temp;
  clear temp
end
