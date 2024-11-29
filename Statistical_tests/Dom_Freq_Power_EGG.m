%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PS : The output matrix has the rows as subjects and columns as the 
% visit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
%addpath('C:\Users\sanke\Downloads\PD_EEG_EGG\fooof_mat-main\fooof_mat');
visits = [1,3];
%subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002 , 1029, 1050 not being considered
subs_to_run = [1001:1067];
filt_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\filtered_EGG_tachy_range\';
PD_Healthy = [repelem(1,34),repelem(0,11),1,1,0,0,0,1,0,0,0,0,0,1,repelem(0,7)]; %1-> patient, 0-> healthy
PD_Healthy([21,24,35,41]) = [];
dom_freq_all = [];
dom_power_all = [];

%% Running the loop to get the maximum power peak

%for electrode = 1:4
  for subs = subs_to_run
     for visit = visits      
       

        %Check if the file exists or not
        if ~isfile(append(filt_path,'V' , num2str(visit),'\',num2str(subs),'.mat'))
            fprintf('%d .mat doesnt exist\n',subs)
            dom_freq(subs-1000,visit) = nan;
            continue
        end
        
        load(append(filt_path,'V' , num2str(visit),'/',num2str(subs),'.mat')); % load the filtered EGG mat file
        
        %{
        % Calculate power
        fs_egg = 250;
        nfft = 60*fs_egg;
        window = hanning(nfft); %rectwin()
        noverlap = nfft/2;
        [Pxx1,w1] =  pwelch(s1_filt,window,noverlap,nfft,fs_egg);
        power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD
        %}

        % This is the PSD by padding the windows for preprandial recording  
        fs_egg = 250;
        windows = floor(size(s1_filt,1) / (60*fs_egg));
        pxx = [];
        for i = 1:windows
            egg_sig = s1_filt(((60*fs_egg)*(i-1))+1:60*fs_egg*i,:); %Taking a window of the signal
            egg_sig_padded = [zeros(2000,4) ; egg_sig ; zeros(2000,4)]; %Zero padding the signal
            %Calculating the psd
            fft_length = size(egg_sig_padded,1);
            [pxx(:,:,i),freq] = pwelch(egg_sig_padded,[],[],fft_length,fs_egg);
            %figure;plot(freq(1:250),pxx_pre(1:250,i));
        end
        psd = mean(pxx,3);
        norm_psd = psd./sum(psd);

        %range_idx = find(freq>0.03 & freq<0.07);
        range_idx = find(freq>0.07 & freq<0.15);
        %range_idx = find(freq>0.0083 & freq<0.03);
        mean_power = mean(norm_psd(range_idx,:),1);
        egg_electrode = find(mean_power == max(mean_power));
        norm_psd = norm_psd(:,egg_electrode);
        psd = psd(:,egg_electrode);

        dom_freq(subs-1000,visit) = freq(norm_psd == max(norm_psd)); 
        dom_power_norm(subs-1000,visit) = max(norm_psd); 
        dom_power(subs-1000,visit) = max(psd);
       

     end    
     
  end

%{
dom_freq([2,29,50],:)= [];
dom_freq([21,24,35,41],:) = [];        
dom_freq(:,4) = mean(dom_freq(:,[1,3]),2);
dom_freq_all(:,electrode) = dom_freq(:,4); %Dominant frequency for all electrodes

dom_power([2,29,50],:)= [];
dom_power([21,24,35,41],:) = [];        
dom_power(:,4) = mean(dom_power(:,[1,3]),2);
dom_power_all(:,electrode) = dom_power(:,4); %Dominant power for all electrodes

end
%}

dom_power_post_minus_pre = dom_power(:,3) - dom_power(:,1);
dom_power_post_minus_pre([2,29,50],:)= [];
dom_power_post_minus_pre([21,24,35,41],:) = []; 
dom_power_post_minus_pre(:,2) = PD_Healthy;
PD = dom_power_post_minus_pre(dom_power_post_minus_pre(:,2)==1);
Healthy = dom_power_post_minus_pre(dom_power_post_minus_pre(:,2)==0);
fprintf('Dominant post-pre power in PD : %f \n',median(PD));
fprintf('Dominant post-pre power in Healthy : %f \n',median(Healthy));
%Wilcoxon test for checking the differences
[p h stats] = ranksum(PD,Healthy)
%{
mean_dom_freq_all = mean(dom_freq_all,2);
mean_dom_freq_all(:,2) = PD_Healthy;
PD = mean_dom_freq_all(mean_dom_freq_all(:,2)==1);
Healthy = mean_dom_freq_all(mean_dom_freq_all(:,2)==0);
fprintf('Dominant frequency for PD (in Hz) : %f \n',mean(PD));
fprintf('Dominant frequency for Healthy (in Hz) : %f \n',mean(Healthy));
%Wilcoxon test for checking the differences
[p h stats] = ranksum(PD,Healthy)


mean_dom_power_all = mean(dom_power_all,2);
mean_dom_power_all(:,2) = PD_Healthy;
PD = mean_dom_power_all(mean_dom_power_all(:,2)==1);
Healthy = mean_dom_power_all(mean_dom_power_all(:,2)==0);
fprintf('Dominant power for PD (in uV^2) : %f \n',mean(PD));
fprintf('Dominant power for Healthy (in uV^2) : %f \n',mean(Healthy));
%Wilcoxon test for checking the differences
[p h stats] = ranksum(PD,Healthy)
%}