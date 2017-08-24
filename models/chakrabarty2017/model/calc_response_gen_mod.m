function [resp,resp_comp]=calc_response_gen_mod(sequence,model,locs)

batchdata=sequence;
numcases = size(batchdata,1);

% if(flag==1)
% %Normalize the data
% data_mean = mean(batchdata,1);
% data_std = std(batchdata);
%
% batchdata =( batchdata - repmat(data_mean,numcases,1) ) ./ ...
%   repmat( data_std, numcases,1);
% end

numhid = model.numhid;     % number of hidden units
numcomp= model.numcomp;      % number of discrete mixture components
numdims = size(batchdata,2); %data (visible) dimension

nt = model.nt;
pastvis=model.pastvis;
pasthid=model.pasthid;
visbiases=model.visbiases;
hidbiases=model.hidbiases;
vishid= model.vishid;
%  pastvis=model.A;
%  pasthid=model.B;
%  visbiases=model.bi;
%  hidbiases=model.bj;
%  vishid= model.w;
min_len=inf;
batchdataindex=nt+1:size(batchdata,1);
nc=length(batchdataindex);
data = single(batchdata(batchdataindex,:));
%data=batchdata;
past = zeros(nc,nt*numdims,'single');

if nargin < 3
  locs = 1:(nt*numdims);
end

count=2;

for hh=nt:-1:1 %note reverse order
               % if(hh>1)
               %   past(:,numdims*(nt-hh)+1:numdims*(nt-hh+1)) =(count*0.3).* batchdata(batchdataindex-hh,:);
               %  count=count-1;
               % else
  past(:,numdims*(nt-hh)+1:numdims*(nt-hh+1)) = batchdata(batchdataindex-hh,:);
  %end

end
% weight_vec=[0.1*ones(1,400) 0.3*ones(1,400) 0.6*ones(1,400)];
% for n=1:nc
%     temp=past(n,:);
%
%     temp=temp.*weight_vec;
%     past(n,:)=temp;
%     clear temp
% end
%Note that we will re-use the effective visible, hidden biases several
%times so we compute them here (per-component) and keep them around
effvisbiases = zeros(nc,numdims,numcomp,'single');
effhidbiases = zeros(nc,numhid,numcomp,'single');
resp=zeros(nc,numhid);

for cc=1:numcomp
  bistar = past*pastvis(locs,locs,cc);
  bjstar = past*pasthid(locs,:,cc);



  effvisbiases(:,:,cc) = repmat(visbiases(cc,locs),nc,1)+ bistar;
  effhidbiases(:,:,cc) = repmat(hidbiases(cc,:),nc,1)+ bjstar;

  resp=resp+(data*vishid(locs,:,cc) + effhidbiases(:,:,cc));
  resp_comp{cc}=data*vishid(locs,:,cc) + effhidbiases(:,:,cc);
  % resp_comp{cc}=data*vishid(:,:,cc) + effhidbiases(:,:,cc);
  %   numcases=size(resp{cc},1);
  %   data_mean1 = mean(resp{cc},1);
  %   data_std1 = std(resp{cc});
  %   resp{cc} =( resp{cc} - repmat(data_mean1,numcases,1) ) ./ ...
  %   repmat( data_std1, numcases,1);
  %   clear data_mean1 data_std1


end
