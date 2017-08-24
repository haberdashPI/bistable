function [resp_data_without_mom]= coherence_wts_comb(data,flag)

w_mod=zeros(size(data,2),size(data,2),size(data,1)+1);
w_mod=weights_hebb_calc(data,w_mod);

resp_data=[];

for j=1:size(w_mod,3)-1
  temp_1=data(j,:)';



  temp_11=w_mod(:,:,j+1)*temp_1;

  if(max(temp_11)>0)
    temp_22=(temp_11./max(temp_11))*2;
  else
    temp_22=-(temp_11./max(temp_11))*2;
  end
  resp_data=[resp_data;temp_22'];
  clear temp_1 temp_11;

end

resp_data_without_mom = resp_data;
