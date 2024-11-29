%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the number of data points of EEG for each task from the mat_files
% and then compare with time stamps on the file names. To segment the EGG
% data, use the event component from the set files. 
% PS : The output cell has the rows as tasks and the column is the
% subject. Inside each cell, the first element is frontal, 2nd element is
% temporo-occipital and 3rd element is centro-parietal. The third dimension
% is the visit (1 and 3).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\sanke\Downloads\Sanket\eeglab2023.1'));
egg_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\EGG_taskwise\';
mat_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\mat_files\';
set_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\';
%filt_path = '/data3/Sanket/filtered_EGG_brady_range/';

subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; % subs 1002,1029 and 1050 not considered
%subs_to_run = [1001,1003:1021,1023:1028,1029:1032,1034:1042,1044:1049,1051:1067]; %sub 1022,1033,1043 not considered due to some issues
%subs_to_run = [1062];
visits = [1,3];
problematic_EGG = []; %The first column is for the subject id, 2nd col -> visit, 3rd col-> task
missing_EGG = []; %1st col -> Visit, 2nd col -> Subject ID 
zero_padded_EGG = []; 

%% Running the loop to get PAC values

for visit = visits

    for subs = subs_to_run

        fileList_egg = dir(append(egg_path,'V',num2str(visit),'/',num2str(subs),'/*.mat'));
        fileList_eeg =  dir(append(set_path,'V',num2str(visit),'/',num2str(subs),'/*.set'));

        for files = 1:4
            
            %Check whether file exists or not
            if ~isfile(append(egg_path,'V' , num2str(visit),'/',num2str(subs),'/Task',num2str(files),'.mat'))
                fprintf('%d .mat doesnt exist\n',subs)
                missing_EGG(end+1,1) = visit;
                missing_EGG(end+1,2) = subs;
                continue
            end

            load(append(egg_path,'V' , num2str(visit),'/',num2str(subs),'/Task',num2str(files),'.mat'));            
            fprintf('EGG file loaded\n')        
       
            %% Finding the EGG electrode with the highest power 
            
            % Filter signal
            fprintf('Being filtered\n')
            s1_filt = bandpass(egg_section,[0.0083 0.03],250);
            fprintf('FIltering done\n')  
            
            %Incase the signal is less than 15000
            if size(s1_filt,1)<15000
                end_point = size(s1_filt,1);
                s1_filt(size(s1_filt,1):15000,:) = 0;
            end
        
            % Calculate power
            fs_egg = 250;
            nfft = 60*fs_egg;
            window = hanning(nfft); %rectwin()
            noverlap = nfft/2;
            [Pxx1,w1] =  pwelch(s1_filt,window,noverlap,nfft,fs_egg);
            power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

            %range_idx = find(w1>0.03 & w1<0.07); %Normogastric range
            %range_idx = find(w1>0.07 & w1<0.15); %Tachygastric range
            range_idx = find(w1>0.0083 & w1<0.03); %Bradygastric range
            %mean_power = mean(power_norm(range_idx,:),1);
            rel_power = sum(power_norm(range_idx,:))./sum(power_norm); %Normogastric range
            egg_electrode = find(rel_power == max(rel_power));

            %Just in case you get two electrodes with the same values :
            %Randomly choose an electrode
            if size(egg_electrode,2)>1
                chosen_idx = randperm(size(egg_electrode,2),1);
                egg_electrode = egg_electrode(chosen_idx);
            end

            if size(s1_filt,1)==15000
                s1_filt(end_point+1:15000,:)= [];
            end
        
            egg_section = s1_filt(:,egg_electrode); %selecting the EGG electrode with the highest relative power
       

            %% get the preprocessed eeg data 

            EEG_data = pop_loadset(append(set_path,'V',num2str(visit),'/',num2str(subs),'/',fileList_eeg(files).name));  
                  
            %% Calculating PAC between beta EEG and the filtered EGG

            task = str2num(fileList_eeg(files,1).name(2)); %get the task number from the filename

            EEG_beta = bandpass(EEG_data.data',[13,30],250); %filtering the EEG signal in beta freq range
            EEG_beta_norm = (EEG_beta-mean(EEG_beta,1))./std(EEG_beta,1); %Normalizing the EEG data (v3 addition)
            EGG_norm = (egg_section-mean(egg_section,"all"))/std(egg_section,[],"all"); %Normalizing the EGG data 

            % Calculating PAC brain region wise

            for region = 1:19
                % get coupled signal
                pac_amp = abs(hilbert(EEG_beta_norm(:,region)));
                egg_phase = angle(hilbert(EGG_norm));
                egg_eeg_coupled_signal = pac_amp.*exp(1i*egg_phase);
                % get top 5 percentile indices of amplitude
                percentile_threshold = prctile(pac_amp, 95);
                indices_top_5_percentile = find(pac_amp >= percentile_threshold);
                % compute overall PAC
                pac_numerator = abs(sum(egg_eeg_coupled_signal));
                pac_denominator = sqrt(sum(pac_amp.*pac_amp))*sqrt(size(EEG_beta_norm,1));
                pac_beta{task,subs-1000,visit}(region) = pac_numerator/pac_denominator;
                % compute top 5 percentile PAC
                pac_numerator = abs(sum(egg_eeg_coupled_signal(indices_top_5_percentile)));
                pac_denominator = sqrt(sum(pac_amp(indices_top_5_percentile).*pac_amp(indices_top_5_percentile)))*sqrt(size(EEG_beta_norm(indices_top_5_percentile),1));
                pac_beta_top5{task,subs-1000,visit}(region) = pac_numerator/pac_denominator;
            end



        end
    end
end

save('beta_PAC_brady_range_19ch_v2.mat','pac_beta','pac_beta_top5'); %Save the variables
