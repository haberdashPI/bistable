clc;
clear all
%addpath('deepmat-master');
addpath('../stimulus');
addpath('../model');

loadload;

paras = [10 8 -2 -1];

base_dir = '/Volumes/Data/Little_Bistable_2017_08_15/deb';
model=load([base_dir '/model.mat']);

generic_file_2 = [base_dir '/generic_mod_newest/'];
%model=load('/Users/debmalyachakrabarty/research/icassp_2016/debcode_maincomp/models/16_06/layer1.mat');
%model=load('/Users/debmalyachakrabarty/research/icassp_2016/debcode_maincomp/models/16_06/model.mat');
%generic_file_2='/Users/debmalyachakrabarty/research/icassp_2016/generic_mod_ekdom_newest/';
%generic_file_2='/Users/debmalyachakrabarty/research/icassp_2016/generic_mod_mixture/';
%locs=load('/Users/debmalyachakrabarty/research/icassp_2016/meet_results/13_07_16/harmonics/comb.mat');
%locs=load('/Users/debmalyachakrabarty/research/icassp_2016/meet_results/13_07_16/onset/onset_stops.mat');

%locs=load('/Users/debmalyachakrabarty/research/icassp_2016/meet_results/13_07_16/models_no/faster_L1.mat');
%baal_neurons=locs.all_neurons;


fs=8000;
tt = [1:0.06*fs]/fs;

% seq_dur=[1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000];



sil=zeros(1,60*8);

deltaf=[1 3 6 9 15];
% sil_1=zeros(1,100*8);

freq_a=[500 750 1000 1250 1500 1750 2000] ;
altered_a=800;

altered_tone=sin(2*pi*800*tt);
% freq_b1=freq_a.*(2.^(deltaf/1));
% freq_b2=freq_a.*(2.^(deltaf/3));
% freq_b3=freq_a.*(2.^(deltaf/6));
% freq_b4=freq_a.*(2.^(deltaf/9));


%tone_a=create_tone(ref_tone,sil);

%tone_orig=create_tone(ref_tone,sil,floor(seq_dur(1)/60));

kuch_bhi_var=1:350;

%dur_sil=[0 2 5 8 11 14 17 20 23 26 29 30];
dur_sil=[0 5 10 20 30 40 50 60 70 80 90 100];
%dur_sil=[0 5 10 15 20 25 30 35 40 45 50 60];
%dur_sil=[0 20 40 60 80 100 120 140 160 180 200];

rand_phase=rand(1,128);

%loc=setdiff(kuch_bhi_var,baal_neurons);
loc=1:350;
% sec_layer=load('/Users/debmalyachakrabarty/research/icassp_2016/meet_results/27_10_16/seq_neurons.mat');

for iii=1:7

  % for jj=1:11
  %         sil1=zeros(1,dur_sil(jj)*8);
  %         sil2=zeros(1,(100-dur_sil(jj))*8);
  %         sil_temp=zeros(1,dur_sil(1)*8);

  for j=1:5
    ref_tone=sin(2*pi*freq_a(iii)*tt);
    freq_b=freq_a(iii).*(2.^(deltaf(j)/12));
    comp_tone=sin(2*pi*freq_b*tt);
    % ref_tone=ripple_create(100,freq_a(iii),freq_a(iii)+5,25,100,rand_phase);
    % comp_tone=ripple_create(100,freq_b,freq_b+5,25,100,rand_phase);
    ref_tone_seq=create_tone(ref_tone,sil,ceil(2000/120));
    comp_tone_seq=create_tone(sil,comp_tone,ceil(2000/120));
    %actual_tone=create_tone(ref_tone,comp_tone,sil,ceil(1000/120));
    % actual_tone_a=create_tone(ref_tone,sil_temp,sil,ceil(1000/240));
    %actual_tone_b=create_tone(comp_tone,sil1,sil2,ceil(1000/240));
    actual_tone=ref_tone_seq+comp_tone_seq;
    actual_tone_alter=[ ref_tone sil  ref_tone sil actual_tone comp_tone sil comp_tone sil];
    actual_tone_orig=[  ref_tone sil  ref_tone sil actual_tone ref_tone sil ref_tone sil];
    baal_tone=[ ref_tone sil ref_tone sil actual_tone];
    y_tone_orig=wav2aud(baal_tone,paras);
    y_tone_actual_orig=wav2aud(actual_tone_orig,paras);
    y_tone_actual_alter=wav2aud(actual_tone_alter,paras);
    layer1_orig{1}=hidden_act(y_tone_orig,model.model);
    %    layer1_orig_avg{1}=zero_mean_unit_var(layer1_orig{1});
    mom_orig{1}=inertia_delta(layer1_orig{1},4);


    %mom_orig{1}= layer1_orig{1};
    mod_orig{1}=mom_orig{1}(:,loc);
    layer1_actual_orig{1}=hidden_act(y_tone_actual_orig,model.model);
    mom_actual_orig{1}=inertia_delta(layer1_actual_orig{1},4);
    %   mom_actual_orig{1}=layer1_actual_orig{1};
    mod_actual_orig{1}=mom_actual_orig{1}(:,loc);
    layer1_actual_alter{1}=hidden_act(y_tone_actual_alter,model.model);
    mom_actual_alter{1}=inertia_delta(layer1_actual_alter{1},4);

    %     mom_actual_alter{1}=layer1_actual_alter{1};
    mod_actual_alter{1}=mom_actual_alter{1}(:,loc);
    for ii=2 %1:8
      generic_model_filename2=[generic_file_2 'mod_ep' int2str(ii) '.mat'];
      generic_mod2=load(generic_model_filename2);
      data_orig{1}=rate_filter( mod_orig{1},ii,size( mod_orig{1},2));
      data_actual_orig{1}=rate_filter(mod_actual_orig{1},ii,size( mod_actual_orig{1},2));
      data_actual_alter{1}=rate_filter(mod_actual_alter{1},ii,size(mod_actual_alter{1},2));
      actual_locs_2nd_layer=loc_find(loc,ii);
      response_orig=calc_response_gen_mod( data_orig{1},generic_mod2,actual_locs_2nd_layer);
      response_orig=avg_layer(response_orig);
      response_actual_orig=calc_response_gen_mod(  data_actual_orig{1},generic_mod2,actual_locs_2nd_layer);
      response_actual_orig=avg_layer(response_actual_orig);
      response_actual_alter=calc_response_gen_mod( data_actual_alter{1},generic_mod2,actual_locs_2nd_layer);
      response_actual_alter=avg_layer(response_actual_alter);
      hebb_orig=coherence_wts_comb(response_orig,0);
      %          %  hebb_orig=hebb_orig(:,sec_layer.select_neurons);
      hebb_actual_orig=coherence_wts_comb(response_actual_orig,0);
      %           % hebb_actual_orig=hebb_actual_orig(:,sec_layer.select_neurons);
      hebb_actual_alter=coherence_wts_comb(response_actual_alter,0);
      % hebb_actual_alter=hebb_actual_alter(:,sec_layer.select_neurons);
      dist_orig_orig(j,ii)=sqrt(trace((hebb_actual_orig(end-5:end,:)-hebb_orig(end-5:end,:))*(hebb_actual_orig(end-5:end,:)-hebb_orig(end-5:end,:))'))
      dist_orig_alter(j,ii)=sqrt(trace((hebb_actual_alter(end-5:end,:)-hebb_orig(end-5:end,:))*(hebb_actual_alter(end-5:end,:)-hebb_orig(end-5:end,:))'))
      %   dist(j,ii)=sqrt(trace((hebb_actual_alter(end-5:end,:)-hebb_actual_orig(end-5:end,:))*(hebb_actual_alter(end-5:end,:)-hebb_actual_orig(end-5:end,:))'))
      %  dist_full_orig(j,ii)=sqrt(trace((hebb_actual_orig-hebb_actual_orig)*(hebb_actual_alter-hebb_actual_orig)'))

      %        dist(j,i,ii)=sqrt(trace((response_actual_alter(end-5:end,:)-response_actual_orig(end-5:end,:))*(response_actual_alter(end-5:end,:)-response_actual_orig(end-5:end,:))'));
      %         dist_full(j,i,ii)=sqrt(trace((response_actual_alter-response_actual_orig)*(response_actual_alter-response_actual_orig)'));

      %         dist_orig_orig(j,ii)=sqrt(trace((response_actual_orig(end-5:end,:)-response_orig(end-5:end,:))*(response_actual_orig(end-5:end,:)-response_orig(end-5:end,:))'))
      %        dist_orig_alter(j,ii)=sqrt(trace((response_actual_alter(end-5:end,:)-response_orig(end-5:end,:))*(response_actual_alter(end-5:end,:)-response_orig(end-5:end,:))'))


      clear generic_model_filename2 generic_mod2 data_actual_orig data_orig response_orig response_actual_alter response_actual_orig hebb_orig hebb_actual_orig hebb_actual_alter
    end

    %  end

  end
  dist_freq_orig(:,:,iii)=dist_orig_orig;
  dist_freq_alter(:,:,iii)=dist_orig_alter;
  iii
end


mean_dist_freq_orig = mean(dist_freq_orig,3);
mean_dist_freq_alter = mean(dist_freq_alter,3);
% % % % % % % % % % %
% % % % % % % %
correct=zeros(5,11);
false=zeros(5,11);
sum_baal_diff_orig=mean_dist_freq_orig(:,2)
sum_baal_diff_alter=mean_dist_freq_alter(:,2)
%   diff_baal=sort(sum_baal_diff_alter-sum_baal_diff_orig)
for ii=1:11
  for n=1:5

    for i=1:100
      if( sum_baal_diff_orig(n)<0.8+0.2*randn(1))
        correct(n,ii)=correct(n,ii)+1;
      end
    end

  end
end

for ii=1:11
  for n=1:5

    for i=1:100
      if( sum_baal_diff_alter(n)<0.8+0.2*randn(1))
        false(n,ii)=false(n,ii)+1;
      end
    end

  end
end
for n=1:5
  for nn=1:11
    [dp_mbd_harm(n,nn),c_mbd_harm(n,nn)] = dprime_simple(correct(n,nn)/101,false(n,nn)/101);
  end
end

mean_hebb_sync1=mean(abs(dp_mbd_harm),2)

%





% % % %
% %  for i=1:5
% %      std_hebb_sync(i)=std(dp_mbd_harm(i,:));
% %  end
% % % %
% % % for i=1:4
% % %     count_sort(i,:)=sort(count(i,:));
% % % end
