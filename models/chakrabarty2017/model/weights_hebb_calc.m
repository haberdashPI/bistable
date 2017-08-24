function w=weights_hebb_calc(x,w)

% corr_grp=zeros(size(x,2),size(x,1),size(x,2));
% for i=1:size(x,2)
%     temp=x(:,i);
%    % temp=temp./max(temp);
%     for j=1:size(x,2)
%         if(i~=j)
%             temp1=x(:,j);
%           %  temp1=temp1./max(temp1);
%             corr_grp(j,:,i)=temp.*temp1;
%
%             clear temp1
%         end
%     end
%
%     clear  temp
%
% end
%
%    for i=1:size(x,2)
%           %w(i,:,:)= HebbNN_temp(squeeze(corr_grp(:,:,i)),squeeze(w(i,:,:)),i);
%           w(i,:,:)= HebbNN_temp(x,squeeze(w(i,:,:)),i);
%           %w_mod(i+1,:)=w_mod(i,:)+w;
%
%    end

for i=1:size(x,1)
  temp=x(i,:);
  if(max(temp)>0)
    temp=temp./max(temp);
  else
    temp=-(temp./max(temp));
  end
  for ii=1:300
    temp_1=temp(ii);
    for j=1:300
      if(j~=ii)
        temp_2=temp(j);
        if(temp_1>=0.9 && temp_2>=0.9)
          %                 if(i<=hist)
          %                 w(ii,j,i+1)= w(ii,j,i)+(temp_1*temp_2*5);
          %                 w(j,ii,i+1)= w(j,ii,i)+temp_1*temp_2*5;
          %                 else

          w(ii,j,i+1)= w(ii,j,i)+temp_1*temp_2;
          w(j,ii,i+1)= w(j,ii,i)+temp_1*temp_2;
          %        end
        elseif(temp_1>0.9 && temp_2<0.9)
          %                 if(i<=hist)
          %                 w(ii,j,i+1)= w(ii,j,i)-temp_1*temp_2*5;
          %                 w(j,ii,i+1)= w(j,ii,i)-temp_1*temp_2*5;
          %                 else
          w(ii,j,i+1)= w(ii,j,i)-temp_1*temp_2;
          w(j,ii,i+1)= w(j,ii,i)-temp_1*temp_2;
          %   end
        elseif(temp_1<0.9 && temp_2>0.9)
          %                 if(i<=hist)
          %                 w(ii,j,i+1)= w(ii,j,i)-temp_1*temp_2*5;
          %                 w(j,ii,i+1)= w(j,ii,i)-temp_1*temp_2*5;
          %                 else
          w(ii,j,i+1)= w(ii,j,i)-temp_1*temp_2;
          w(j,ii,i+1)= w(j,ii,i)-temp_1*temp_2;
          %    end
          %             elseif(temp_1<=0 || temp_2>0)
          % %                  w(j,ii,i+1)= w(j,ii,i)-0.1;
        else
          %                  if(i<=hist)
          %                 w(ii,j,i+1)= w(ii,j,i)+temp_1*temp_2*5;
          %                 w(j,ii,i+1)= w(j,ii,i)+temp_1*temp_2*5;
          %                 else

          w(ii,j,i+1)= w(ii,j,i)+temp_1*temp_2;
          w(j,ii,i+1)= w(j,ii,i)+temp_1*temp_2;
          %    end
        end
      end
    end
  end
end

% for i=1:400
%     temp=squeeze(w(i,:,:));
%     for ii=1:400
%         temp1=temp(ii,:);
%         if(max(temp1)~=0)
%         temp1=temp1./max(temp1);
%         end
%         w_mod(i,ii,:)=temp1;
%         clear temp1
%     end
%     clear temp
% end
