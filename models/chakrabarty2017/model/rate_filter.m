      
function data=rate_filter(data_act,hist,dims)







%          if(hist==1)
%              data=data_act;
%              
%          else
%              count=1;
%              for t=1:hist:size(data_act,1)-hist
%                  temp=data_act(t:t+hist-1,:);
%                  temp1=temp(:);
%                  data(count,:)=temp1;
%                  clear temp temp1
%                  count=count+1;
%              end
% 
%              
%          end

 if(hist==1)
             data=data_act;
             
 else
         temp=data_act';
         temp=temp(:,1:floor(size(temp,2)/hist)*hist);
         temp1=reshape(temp,dims*hist,floor(size(temp,2)/hist));
         data=temp1';

             
         end