clc;
clear all
%addpath('/Users/debmalyachakrabarty/research/icassp_2016/debcode_maincomp/code/');
addpath('../stimulus');
addpath('../model');


loadload;

paras = [4 4 -2 -1];

base_dir = '/Volumes/Data/Little_Bistable_2017_08_15/deb';
model=load([base_dir '/model.mat']);

generic_file_2 = [base_dir '/generic_mod_newest/'];

 fs=8000;
 tt = [1:0.125*fs]/fs;

 seq_dur=[500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000];

 sil=zeros(1,125*8);

 deltaf=[1 3 6 9];
% sil_1=zeros(1,100*8);

freq_a=[500 1000 2000 4000] ;
altered_a=800;

%altered_tone=sin(2*pi*800*tt);
% freq_b1=freq_a.*(2.^(deltaf/1));
% freq_b2=freq_a.*(2.^(deltaf/3));
% freq_b3=freq_a.*(2.^(deltaf/6));
% freq_b4=freq_a.*(2.^(deltaf/9));


%tone_a=create_tone(ref_tone,sil);

%tone_orig=create_tone(ref_tone,sil,floor(seq_dur(1)/60));
% sec_layer=load('/Users/debmalyachakrabarty/research/icassp_2016/meet_results/27_10_16/seq_neurons.mat');

for iii=1:1
for i=1:19
    for j=1:4
        ref_tone=sin(2*pi*freq_a(iii)*tt);
        freq_b=freq_a(iii).*(2.^(deltaf(j)/12));
        comp_tone=sin(2*pi*freq_b*tt);
     %   ref_tone_seq=create_tone(ref_tone,sil,floor(seq_dur(i)/60));
      %  comp_tone_seq=create_tone(sil,comp_tone,floor(seq_dur(i)/60));

        orig_tone=[ref_tone comp_tone ref_tone];

        actual_tone=create_tone(orig_tone,sil,ceil(seq_dur(i)/500));
      %  actual_tone=ref_tone_seq+comp_tone_seq;
        actual_tone_alter=[ref_tone sil ref_tone sil actual_tone comp_tone sil comp_tone sil];
        actual_tone_orig=[ref_tone sil ref_tone sil actual_tone ref_tone sil ref_tone sil];
      %  y_tone_orig=wav2aud(tone_orig,paras);
        y_tone_actual_orig=wav2aud(actual_tone_orig,paras);
        y_tone_actual_alter=wav2aud(actual_tone_alter,paras);
      %  layer1_orig{1}=hidden_act(y_tone_orig,model.model);
       % mom_orig{1}=inertia_delta(layer1_orig{1},4);
        layer1_actual_orig{1}=hidden_act(y_tone_actual_orig,model.model);
         mom_actual_orig{1}=inertia_delta(layer1_actual_orig{1},4);
        layer1_actual_alter{1}=hidden_act(y_tone_actual_alter,model.model);
         mom_actual_alter{1}=inertia_delta(layer1_actual_alter{1},4);
        for ii=1:8
           generic_model_filename2=[generic_file_2 'mod_ep' int2str(ii) '.mat'];
           generic_mod2=load(generic_model_filename2);
          % data_orig{1}=rate_filter( layer1_orig{1},ii,350);
           data_actual_orig{1}=rate_filter( layer1_actual_orig{1},ii,350);
           data_actual_alter{1}=rate_filter(layer1_actual_alter{1},ii,350);
          % response_orig=calc_response_gen( data_actual_orig{1},generic_mod2);
           response_actual_orig=calc_response_gen_mod(  data_actual_orig{1},generic_mod2);
           response_actual_alter=calc_response_gen_mod( data_actual_alter{1},generic_mod2);
         %  hebb_orig=coherence_wts_comb(response_orig,1);
           hebb_actual_orig=coherence_wts_comb(response_actual_orig,1);%hebb_actual_orig=hebb_actual_orig(:,sec_layer.select_neurons);
           hebb_actual_alter=coherence_wts_comb(response_actual_alter,1);%hebb_actual_alter=hebb_actual_alter(:,sec_layer.select_neurons);
           %dist_orig_orig(j,i,ii)=sqrt(trace((hebb_actual_orig(end-5:end,:)-hebb_orig(end-5:end,:))*(hebb_actual_orig(end-5:end,:)-hebb_orig(end-5:end,:))'));
           %dist_orig_alter(j,i,ii)=sqrt(trace((hebb_actual_alter(end-5:end,:)-hebb_orig(end-5:end,:))*(hebb_actual_alter(end-5:end,:)-hebb_orig(end-5:end,:))'));
           dist(j,i,ii)=sqrt(trace((hebb_actual_alter(end-5:end,:)-hebb_actual_orig(end-5:end,:))*...
                                   (hebb_actual_alter(end-5:end,:)-hebb_actual_orig(end-5:end,:))'))
          % dist_full(j,i,ii)=sqrt(trace((hebb_actual_alter-hebb_actual_orig)*(hebb_actual_alter-hebb_actual_orig)'));

         %  dist(j,i,ii)=sqrt(trace((response_actual_alter(end-5:end,:)-response_actual_orig(end-5:end,:))*(response_actual_alter(end-5:end,:)-response_actual_orig(end-5:end,:))'));
         %  dist_full(j,i,ii)=sqrt(trace((response_actual_alter-response_actual_orig)*(response_actual_alter-response_actual_orig)'));

           clear generic_model_filename2 generic_mod2 data_actual_orig data_orig response_orig response_actual_alter response_actual_orig hebb_orig hebb_actual_orig hebb_actual_alter
        end

    end
   i
end

%sum_dist=squeeze(sum(dist,3));
dist_freq(:,:,:,iii)=dist;
%sum_dist_full=squeeze(sum(dist_full,3));
%dist_freq_full(:,:,iii)=sum_dist_full;
clear dist sum_dist dist_full sum_dist_full
iii
end

% % avg_dist=squeeze(mean(dist_freq,4));
  avg_dist_main = sum(dist_freq(:,:,1:3),3);
 %%% small trick to interprete the results in terms of buildup
 for i=1:4
     temp(i,:)=avg_dist_main(i,1)-avg_dist_main(i,:);
 end
% % % % % % %
% % % %  %  mean_dist_freq=mean(dist_freq_full,3);
% % % % % % %
% % % %
count=zeros(4,19,100);
for ii=1:100
for n=1:4
    for nn=1:19
       for i=1:100
           if(temp(n,nn)>4+1*randn(1))
               count(n,nn,ii)=count(n,nn,ii)+1;
           end
       end
    end
end
end
mean_count=mean(count,3)./100
% % for i=1:4
% %     count_sort(i,:)=sort(mean_count(i,:));
% % end
