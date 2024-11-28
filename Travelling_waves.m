%% Loading some files
addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
%addpath(genpath('C:\Users\sanke\Downloads\klabhub-bayesFactor-04b80fd'));
filepath = 'C:\Users\sanke\Downloads\PD_EEG_EGG\filtered_EGG_tachy_range\V1\';
filepath_EEG = 'C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\V1\';

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

%% Calculating travelling waves for EGG data

% PS : To change the order of the electrodes go down to the functions

%Getting different cohorts
PD = subject_info(subject_info(:,2)==1,1);
Healthy = subject_info(subject_info(:,2)==0,1);
dysPD = subject_info(subject_info(:,5)>0,1);
nondys = subject_info(subject_info(:,5)==0,1);

% For PD subjects
[FW_PD, BW_PD] = logratio_EGG(PD, filepath);

% For Healthy subjects
[FW_healthy, BW_healthy] = logratio_EGG(Healthy, filepath);

% For dyskinesia subjects
[FW_dys, BW_dys] = logratio_EGG(dysPD, filepath);

% For nondyskinesia subjects
[FW_nondys, BW_nondys] = logratio_EGG(nondys, filepath);

% Calculating Log ratios
logratio_PD = log(FW_PD./BW_PD);
logratio_healthy = log(FW_healthy./BW_healthy);
logratio_nondys = log(FW_nondys ./ BW_nondys);

% One sample wilcoxon test (Hypothesis : Is the median different than zero)
[p_PD, ~] = signrank(logratio_PD);
[p_healthy, ~] = signrank(logratio_healthy);
[p_nondys, ~] = signrank(logratio_nondys);

%% Computing travelling waves for EEG waves

%Getting different cohorts
PD = subject_info(subject_info(:,2)==1,1);
Healthy = subject_info(subject_info(:,2)==0,1);
dysPD = subject_info(subject_info(:,5)>0,1);

% For PD subjects
[FW_PD, BW_PD] = logratio_EEG(PD, filepath_EEG);


% For Healthy subjects
[FW_healthy, BW_healthy] = logratio_EEG(Healthy, filepath_EEG);

% Calculating Log ratios
logratio_PD = log(FW_PD./BW_PD);
logratio_healthy = log(FW_healthy./BW_healthy);

% One sample ttest (Hypothesis : Is the mean different than zero)
[bf_PD, p_PD] = bf.ttest(logratio_PD);
[bf_healthy, p_healthy] = bf.ttest(logratio_healthy);

% Two sample ttest
[h pval] = ttest2(logratio_PD, logratio_healthy);

%% Reserve code for plotting 


load(append(filepath,'\', num2str(1053),'.mat'));
egg_sig(:,1) = s1_filt(:,4)';
egg_sig(:,2) = s1_filt(:,3)';
egg_sig(:,3) = s1_filt(:,2)';

num_epochs = size(egg_sig,1) / (60*250); % Divided by 60s*250Hz (250 Hz is the sampling rate)

% THe loop for creating 
for epoch = 1:num_epochs

egg_sig_1min = egg_sig(1+(epoch-1)*60*250:epoch*60*250,:)';      
%}

a=1; %width pixel to compute the frequency
[m,n]=size(egg_sig_1min); %numberPixel x and y
fx=(1/a)*((1:n)-mean(1:n));
fy=(1/a)*((1:m)-mean(1:m));
twod_fft = abs(fftshift(fft2(egg_sig_1min)));
[thX1,thX2]=findingXboundaries(twod_fft,fx);

figure
subplot(2,1,1)
surf(egg_sig_1min)
shading interp; view(0,90); axis tight; colorbar;
xlabel('time [ms]')
ylabel('layers')
title('Looking for travelling waves')
subplot(2,2,3)
surf(fx,fy,twod_fft)
title('2D FFT - zoomed xAxis')
xlabel('timeFreq')
ylabel('layersFreq')
shading interp; view(0,90); axis tight; colorbar;
set(gca,'xlim',[thX1 thX2],'clim',[min(min(twod_fft)) max(max(twod_fft))],'yscale','lin');
subplot(2,2,4)
twod_fft(twod_fft==0)=NaN;
upperPart=twod_fft(fy>0,(fx>0 & fx<thX2));
lowerPart=twod_fft(fy<0,(fx>0 & fx<thX2));
bar([1 2],[mean(upperPart(:)) mean(lowerPart(:))])
title('Quantifying the waves')
xlabel('upper right quadrant x>0, y>0 -  lower right quadrant x>0, y<0')
ylabel('sum of the quadrant')
sgtitle(epoch)

end

%% Function to find X boundaries

function [thX1,thX2]=findingXboundaries(twod_fft,fx)

    %this is just a makiavelic way to find the boundaries of the interesting frequencies (i.e. thX1,thX2) ##to improve
    sumFFT=sum(twod_fft);
    threshold=0.05*max(sumFFT);
    [aa,var]=sort(abs(sumFFT-threshold));
    var=var(1:40); %40 big enough not to cut. ##to improve
    thX1=fx(min(var));
    thX2=fx(max(var));

end

%% Function to find logratios in EGG data 

function [FW, BW] = logratio_EGG(dysPD,filepath)

   count = 1;
   for subs = dysPD'
       %Load subject data
       load(append(filepath,'\', num2str(subs),'.mat'));
       
       % Getting the electrode data : Forward(first electrode) to
       % backward(last electrode)
       egg_sig = [];
       egg_sig(:,1) = s1_filt(:,4);
       egg_sig(:,2) = s1_filt(:,3);
       egg_sig(:,3) = s1_filt(:,2);

       % Calculating number of 1 min epochs 
       num_epochs = size(egg_sig,1) / (60*250); % Divided by 60s*250Hz (250 Hz is the sampling rate)
    
       tempFW = [];
       tempBW = [];

       for epoch = 1:num_epochs
           
           % Computing 2d fft
           egg_sig_1min = egg_sig(1+(epoch-1)*60*250:epoch*60*250,:)';
           twod_fft = abs(fftshift(fft2(egg_sig_1min)));

           a=1; %width pixel to compute the frequency

           [m,n]=size(egg_sig_1min); %numberPixel x and y

           fx=(1/a)*((1:n)-mean(1:n));
           fy=(1/a)*((1:m)-mean(1:m));

           [thX1,thX2]=findingXboundaries(twod_fft,fx);

           upperPart=twod_fft(fy>0,(fx>0 & fx<thX2));
           lowerPart=twod_fft(fy<0,(fx>0 & fx<thX2));

           tempFW(epoch) = mean(upperPart);
           tempBW(epoch) = mean(lowerPart);

       end

       FW(count) = mean(tempFW);
       BW(count) = mean(tempBW);
       count = count + 1;

   end
    
end


%% Function to find logratios for EEG data

function [FW, BW] = logratio_EEG(cohort,filepath)

   count = 1;   
   for subs = cohort'
       
       %Load subject data
       s = ls(fullfile(append(filepath, '\',num2str(subs),'\'),'T4*'));
    
       if ~isempty(s)
           eeg = pop_loadset(append(filepath,'\', num2str(subs),'\',s));
       else 
           continue
       end

       eeg = pop_loadset(append(filepath,'\', num2str(subs),'\',s));

       eeg_sig = [];
    
       %Getting the channel indices for the occipital to frontal axis
       occ2frontal = {'fp1','fz','cz','pz','o1'};
       %occ2frontal = {'t3','c3','cz','c4','t4'}; %left to right hemisphere
       EEG_preproc = pop_loadset('C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\V1\1001\T1(11-25-07).set'); %Load just one EEG set file

       for i=1:size(occ2frontal,2)
           chan_idx(i) = find(strcmpi({EEG_preproc.chanlocs.labels},occ2frontal{i}));
       end
       
       % Filering in the alpha range
       eeg_sig = bandpass(eeg.data(chan_idx,:)', [8 12],250);
       
       %Calculating number of 1s epochs
       num_epochs = size(eeg_sig,1) / (1*250); % Divided by 1s*250Hz (250 Hz is the sampling rate)
    
       tempFW = [];
       tempBW = [];

       for epoch = 1:num_epochs
           % Computing 2d fft
           eeg_sig_1s = eeg_sig(1+(epoch-1)*1*250:1*epoch*250,:)';
           twod_fft = abs(fftshift(fft2(eeg_sig_1s)));

           a=1; %width pixel to compute the frequency

           [m,n]=size(eeg_sig_1s); %numberPixel x and y

           fx=(1/a)*((1:n)-mean(1:n));
           fy=(1/a)*((1:m)-mean(1:m));

           [thX1,thX2]=findingXboundaries(twod_fft,fx);

           upperPart=twod_fft(fy>0,(fx>0 & fx<thX2));
           lowerPart=twod_fft(fy<0,(fx>0 & fx<thX2));

           tempFW(epoch) = nanmean(upperPart(:));
           tempBW(epoch) = nanmean(lowerPart(:));
    
       end

       FW(count) = nanmean(tempFW);
       BW(count) = nanmean(tempBW);
       count = count + 1;
   end
end