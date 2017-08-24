function [temp]=hidden_act(y,R)
%addpath('/Users/debmalyachakrabarty/research/icassp_2016/deepmat-master/');
count=1;
for i=1:3:size(y,1)-3
          temp=y(i:i+2,:);

          temp_vec=temp';
        %  mixture_spect(count,:)=temp_vec(:)./R.sigmas.^2;
          mixture_spect(count,:)=temp_vec(:);
         count=count+1;
         clear temp_vec temp


end
%
         patches_mean = mean(mixture_spect, 2);
         mixture_spect = bsxfun(@minus, mixture_spect, patches_mean);
         patches_std = std(mixture_spect, [], 2);
         mixture_spect = bsxfun(@rdivide, mixture_spect, patches_std);

         % size(mixture_spect)

%        temp=bsxfun(@plus,bsxfun(@rdivide,mixture_spect,R.sigmas.^2)*R.W, R.hbias');
        temp=bsxfun(@plus,mixture_spect*R.W{1}, R.biases{2}');
     %   h=1./(1+exp(-temp));
      %  res1=binornd(1,h,size(h,1),size(h,2));
%        for j=1:size(h,1)
%            temp=h(j,:);
%            res1(j,:) = double(temp > rand(1,400)); %Activating hiddens
%        end
