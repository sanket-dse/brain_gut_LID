%% Adding some files to the path 

addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
subs_to_run = [1001:1067];
%subs_to_run = [1001,1006];
visits = [1,3];
normo_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\filtered_EGG_normo_range\';
brady_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\filtered_EGG_brady_range\';
tachy_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\filtered_EGG_tachy_range\';
electrode = 3;

psd_brady = nan(9501,length(subs_to_run),3);
psd_normo = nan(9501,length(subs_to_run),3);
psd_tachy = nan(9501,length(subs_to_run),3);

%% Calculating the PSD for different EGG ranges

% Loading up EGG ranges files and plotting them 

for subs = subs_to_run

    for visit = visits
                  
        if isfile(append(brady_path,'V',num2str(visit),'\',num2str(subs),'.mat'))
            load(append(brady_path,'V',num2str(visit),'\',num2str(subs),'.mat'))
            [~,psd_brady(:,subs-1000,visit), freq] = paddedPSD(s1_filt, electrode);
        
        else
            psd_brady(:,subs-1000,visit) = nan(9501,1);
        end

        if isfile(append(normo_path,'V',num2str(visit),'\',num2str(subs),'.mat'))
            load(append(normo_path,'V',num2str(visit),'\',num2str(subs),'.mat'))
            [~,psd_normo(:,subs-1000,visit), freq] = paddedPSD(s1_filt, electrode);
        else
            psd_normo(:,subs-1000,visit) = nan(9501,1);
        end

        if isfile(append(tachy_path,'V',num2str(visit),'\',num2str(subs),'.mat'))
            load(append(tachy_path,'V',num2str(visit),'\',num2str(subs),'.mat'))
            [~,psd_tachy(:,subs-1000,visit), freq] = paddedPSD(s1_filt, electrode);
        else
            psd_tachy(:,subs-1000,visit) = nan(9501,1);
        end
        
    end
end

psd_brady = nanmean(nanmean(psd_brady,2),3);
psd_normo = nanmean(nanmean(psd_normo,2),3);
psd_tachy = nanmean(nanmean(psd_tachy,2),3);

figure;
plot(freq(1:12),psd_brady(1:12),'LineWidth',2);hold on;
plot(freq(1:12),psd_normo(1:12),'LineWidth',2); 
plot(freq(1:12),psd_tachy(1:12),'LineWidth',2); 
xline(freq(psd_normo==max(psd_normo)),'r--','LineWidth',2);
xline(freq(psd_brady==max(psd_brady)),'b--','LineWidth',2);
xline(freq(psd_tachy==max(psd_tachy)),'y--','LineWidth',2); hold off;
xlabel('Frequencies');
ylabel('Power spectral density(uV^2)');
legend('Bradygastric','Normogastric','Tachygastric');
title('Power spectral densities for different ranges of EGG');


%% Reading the UPDRS reports of subjects 

% Reading the excel file
path1 = 'C:\Users\sanke\Downloads\PD_EEG_EGG\NURO Record.xlsx';
[~,sheets] = xlsfinfo(path1);
[~, ~, patient_record] = xlsread(path1,sheets{:,1});
[~, ~, UPDRS] = xlsread(path1,sheets{:,3});

% Extracting data from UPDRS score sheet

addpath('C:\Users\sanke\Downloads\PD_EEG_EGG\cell2float');
UPDRS_clean = cell2float(UPDRS);
all_scores = nan(64,74);  

for i=1:66
    for j = 2:75
        all_scores(i,j-1) = UPDRS_clean(i+4,j);
     end
end
all_scores(:,30:32) = [];

% Removing subs 1002 and 1050 (since they are not considered for analysis)

idx = find(all_scores(:,1)==1002 | all_scores(:,1)==1050);
all_scores = all_scores([1,3:48,50:66],:);
dys_score = all_scores(:,66) + all_scores(:,65);
dys_score(isnan(dys_score)) = 0 ; 

% Getting subject info from the xlsx file
subject_info = nan(64,4); %1st col -> ID, 2nd col -> Patient/Healthy, 3rd col -> Age, 4th col -> Gender
for i=1:66

    subject_info(i,1) = patient_record{i+2,2};
    temp = patient_record{i+2,5};
    
    if strcmpi(temp,'patients')
        subject_info(i,2) = 1;
    else
        subject_info(i,2) = 0;
    end

    subject_info(i,3) = patient_record{i+2,7};
    subject_info(i,4) = patient_record{i+2,8};
end

% NOTE : Here in 4th col, 77 -> M , 70 -> F 
% In 2nd col, 1 -> Patient ,  0 -> Healthy 
% 5th col -> Section 4.1+4.2 -> Dyskinesia Score
% 1st col -> Subject ID
% 3rd col -> Age

% Removing subs 1002 and 1050 (since they are not considered for analysis)

idx = find(subject_info(:,1)==1002 | subject_info(:,1)==1050);
subject_info(idx,:) = nan;
subject_info = rmmissing(subject_info);
subject_info(:,5) = dys_score;

% Removing subjects that didn't have pre/post pac features at all (Subs 1022, 1025, 1037,
% 1043)
idx = find(subject_info(:,1)==1022 | subject_info(:,1)==1025 | subject_info(:,1)==1037 | subject_info(:,1)==1043);
subject_info(idx,:) = nan;
subject_info = rmmissing(subject_info);

%Getting the subject IDs for healthy and PD subjects
PD_subs = subject_info(subject_info(:,2)==1,1)';
healthy_subs = subject_info(subject_info(:,2)==0,1)';

%% Calculating PSD for healthy and PD patients for pre and post conditions

% Loading up post medication normogastric EGG for PD and healthy 

psd_PD_pre = nan(9501,67);
psd_PD_post = nan(9501,67);
psd_healthy_pre = nan(9501,67);
psd_healthy_post = nan(9501,67);

for sub = PD_subs
    if isfile(append(normo_path,'V',num2str(1),'\',num2str(sub),'.mat'))
        load(append(normo_path,'V',num2str(1),'\',num2str(sub),'.mat'))
        [~,psd_PD_pre(:,sub-1000), freq] = paddedPSD(s1_filt, electrode);
    end
end

for sub = PD_subs    
    if isfile(append(normo_path,'V',num2str(3),'\',num2str(sub),'.mat'))
        load(append(normo_path,'V',num2str(3),'\',num2str(sub),'.mat'))
        [~,psd_PD_post(:,sub-1000), freq] = paddedPSD(s1_filt, electrode);    
    end
end

for sub = healthy_subs    
    if isfile(append(normo_path,'V',num2str(1),'\',num2str(sub),'.mat'))
        load(append(normo_path,'V',num2str(1),'\',num2str(sub),'.mat'))
        [~,psd_healthy_pre(:,sub-1000), freq] = paddedPSD(s1_filt, electrode);    
    end
end

for sub = healthy_subs    
    if isfile(append(normo_path,'V',num2str(3),'\',num2str(sub),'.mat'))
        load(append(normo_path,'V',num2str(3),'\',num2str(sub),'.mat'))
        [~,psd_healthy_post(:,sub-1000), freq] = paddedPSD(s1_filt, electrode);    
    end
end
        
% Taking mean of the PSD 

psd_PD_pre = nanmean(psd_PD_pre,2);
psd_PD_post = nanmean(psd_PD_post,2);
psd_healthy_pre = nanmean(psd_healthy_pre,2);
psd_healthy_post = nanmean(psd_healthy_post,2);

% Plotting the PSDs

figure;
plot(freq(1:12),psd_PD_pre(1:12),'LineWidth',2);hold on;
plot(freq(1:12),psd_healthy_pre(1:12),'LineWidth',2);
plot(freq(1:12),psd_PD_post(1:12),'LineWidth',2);
plot(freq(1:12),psd_healthy_post(1:12),'LineWidth',2);

plot(freq(psd_healthy_pre==max(psd_healthy_pre)),psd_healthy_pre(psd_healthy_pre==max(psd_healthy_pre)),'r*');
plot(freq(psd_PD_pre==max(psd_PD_pre)),psd_PD_pre(psd_PD_pre==max(psd_PD_pre)),'r*');
plot(freq(psd_PD_post==max(psd_PD_post)),psd_PD_post(psd_PD_post==max(psd_PD_post)),'r*');
plot(freq(psd_healthy_post==max(psd_healthy_post)),psd_healthy_post(psd_healthy_post==max(psd_healthy_post)),'r*'); hold off

xlabel('Frequencies');
ylabel('Power spectral density (uV^2)');
legend('pre condition PD','pre condition Healthy','post condition PD','post condition Healthy');
title('Power spectral densities of normogastric EGG for PD and healthy subjects');


%% Function to calculate the PSD 

function [norm_psd, psd, freq] = paddedPSD(s1_filt, electrode)

   fs_egg = 250;
   windows = floor(size(s1_filt,1) / (60*fs_egg));
   pxx = [];
   for i = 1:windows
       egg_sig = s1_filt(((60*fs_egg)*(i-1))+1:60*fs_egg*i,electrode); %Taking a window of the signal
       egg_sig_padded = [zeros(2000,1) ; egg_sig ; zeros(2000,1)]; %Zero padding the signal
            
       %egg_sig = s1_filt(((60*fs_egg)*(i-1))+1:60*fs_egg*i,:); %Taking a window of the signal
       %egg_sig_padded = [zeros(2000,4) ; egg_sig ; zeros(2000,4)]; %Zero padding the signal

       %Calculating the psd
       fft_length = size(egg_sig_padded,1);
       [pxx(:,:,i),freq] = pwelch(egg_sig_padded,[],[],fft_length,fs_egg);
       %figure;plot(freq(1:250),pxx_pre(1:250,i));
   end

   psd = mean(pxx,3);
   norm_psd = psd./sum(psd);

end