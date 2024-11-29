%% Loading the files and adding the path to folders

addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
addpath(genpath('C:\Users\sanke\Downloads\PD_EEG_EGG\PACT-master'));
addpath(genpath('C:\Users\sanke\Downloads\PD_EEG_EGG\spider_plot'));

% Healthy subject
load('C:\Users\sanke\Downloads\PD_EEG_EGG\EGG_tachy_taskwise\V1\1066\Task3.mat');
egg_sig_Healthy = egg_section;

eeg_Healthy = pop_loadset('C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\V1\1066\T3(12-36-08).set');

%PD subject
load('C:\Users\sanke\Downloads\PD_EEG_EGG\EGG_tachy_taskwise\V1\1049\Task3.mat');
egg_sig_PD = egg_section;

eeg_PD = pop_loadset('C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\V1\1049\T3(12-21-17).set');

%PD subject with no dyskinesia 
load('C:\Users\sanke\Downloads\PD_EEG_EGG\EGG_normo_taskwise\V1\1054\Task3.mat');
egg_sig_nondys = egg_section;

eeg_nondys = pop_loadset('C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\V1\1054\T3(11-54-48).set');


%% Calculating coupled signal

%Getting the eeg signal and bandpass filtering it
eeg_sig = eeg_PD.data(3,:); %Just T5 electrode
eeg_sig = bandpass(eeg_sig',[13,30],250);

eegH = hilbert(eeg_sig);
eggH = hilbert(egg_sig_PD);

eeg_amp = abs(eegH);
egg_phase = angle(eggH);
eeg_egg_coupled_signal = eeg_amp.*exp(1i*egg_phase);

percentile_threshold = prctile(eeg_amp, 95);
indices_top_5_percentile = find(eeg_amp >= percentile_threshold);

%Top 5%ile PAC
pac_numerator = abs(sum(eeg_egg_coupled_signal(indices_top_5_percentile)));
pac_denominator = sqrt(sum(eeg_amp(indices_top_5_percentile).*eeg_amp(indices_top_5_percentile)))*sqrt(length(eeg_amp(indices_top_5_percentile)));
pac_top5percentile = pac_numerator/pac_denominator;
%pac_top5percentile = pac_numerator / sqrt(sum(eeg_amp(indices_top_5_percentile).*eeg_amp(indices_top_5_percentile)));

new_coupled_signal = eeg_egg_coupled_signal(indices_top_5_percentile);
new_egg_phase = egg_phase(indices_top_5_percentile);
new_eeg_amp = eeg_amp(indices_top_5_percentile);

%% Visualizing PAC (Phase sorted amplitude bins)

% sorting amplitudes in phase bins
nbins = 49;
coupled_amp = abs(new_coupled_signal)/pac_denominator;
phase_bins = linspace(-pi,pi,nbins+1);
phase_sorted_Amp = []; %first column -> mean amp ; second column -> std amp

for i = 1:nbins

    idx = find(new_egg_phase > phase_bins(i) & new_egg_phase < phase_bins(i+1));
    if isempty(idx)
        amp=0;
    else
        amp = coupled_amp(idx);
    end

    phase_sorted_Amp(end+1,1) = sum(amp); %Sum of amplitude
    %phase_sorted_Amp(end,2) = std(amp); %Std amplitude
end

%{
%plot the phase sorted amplitude histogram
figure;
%subplot(1,2,1)
bar(phase_sorted_Amp(:,1), 'FaceColor', [0 86/255 69/255], 'BarWidth', 1); % Billiard 
hold on
%errorbar(1:nbins, phase_sorted_Amp(:,1), zeros(1,nbins), phase_sorted_Amp(:,2), 'LineStyle','none','Color',[0 0 0])
yMax = max(phase_sorted_Amp(:,1) + phase_sorted_Amp(:,2));
yMin = min(phase_sorted_Amp(:,1));
%yMin = min(phase_sorted_Amp(:,1));
ylim([0.65*max(phase_sorted_Amp(:,1)),1.25*max(phase_sorted_Amp(:,1))])

title('PD dyskinesia subject : Pre-medication normogastric PAC for Task 2 (Fp2 channel)')
%}

%plot the phase sorted amplitude histogram
figure;
%subplot(1,2,1)
bar(phase_sorted_Amp, 'FaceColor', [0 86/255 69/255], 'BarWidth', 1); % Billiard 
hold on
%errorbar(1:nbins, phase_sorted_Amp(:,1), zeros(1,nbins), phase_sorted_Amp(:,2), 'LineStyle','none','Color',[0 0 0])
yMax = max(phase_sorted_Amp);
yMin = min(phase_sorted_Amp);
%ylim([0.65*yMax,1.25*yMax])

title('Healthy subject : Pre-medication tachygastric PAC for Task 3 (C4 channel)')

% overplot a cosine wave
t = 0:0.1:nbins;
oneCycleConstant = nbins/6.28; % 0 to 2pi
%plot(t, (cos(t/oneCycleConstant)+1)/4*(yMax*1-yMin*0.3) + 0.03, 'r', 'LineWidth',2)
plot(t, (cos(t/oneCycleConstant)+1)/4*(yMax*2-yMin*1), 'r', 'LineWidth',2)
set(gca, 'XTickLabel',{'0','pi','2*pi'},'XTick',[0 nbins/2 nbins], 'CLim',[1 2])
xlabel('Phase')
ylabel('Normalized Amplitude')
hold off;

%% Visualizing PAC (Polar plot)
subplot(1,2,2)
coupled_phase = angle(new_coupled_signal);
circ_plot(coupled_phase,'hist',[], nbins, true,true); %'linewidth',0.001)
sgtitle('Healthy subject : Pre-medication normogastric PAC for Task 2 (Fp2 channel)');
%}

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

% Sorting the subjects according to their condition
dysPD = subject_info(subject_info(:,5)>0,1)-1000;
nondysPD = subject_info(subject_info(:,5)==0 & subject_info(:,2)==1,1)-1000;
healthy = subject_info(subject_info(:,5)==0 & subject_info(:,2)==0,1)-1000;

%% Spider plot for different subjects
 
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG_2nd_run\beta_PAC_normo_range_19ch_2.mat');
pac_normo = pac_beta_top5;
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG_2nd_run\beta_PAC_tachy_range_19ch_2.mat');
pac_tachy = pac_beta_top5;
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG_2nd_run\beta_PAC_brady_range_19ch_2.mat');
pac_brady = pac_beta_top5;
clear('pac_beta_top5','pac_beta');
heart_features = readtable('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\Demographics_heart_features_v2.csv');

% Getting task averaged pre HRV and pre HFLF
pre_HRV = table2array(mean(heart_features(:,[1,2,3,4]),2));
pre_HFLF = table2array(mean(heart_features(:,[5,6,7,8]),2));

% Adding the removed subjects as elements to the array
pre_HRV_new = [pre_HRV(1)', nan, pre_HRV(2:20)', nan, pre_HRV(21:22)', nan, pre_HRV(23:25)', nan, pre_HRV(26:32)', nan, pre_HRV(33:37)', nan, ...
    pre_HRV(38:43)', nan, pre_HRV(44:60)'];
pre_HFLF_new = [pre_HFLF(1)', nan, pre_HFLF(2:20)', nan, pre_HFLF(21:22)', nan, pre_HFLF(23:25)', nan, pre_HFLF(26:32)', nan, pre_HFLF(33:37)', nan, ...
    pre_HFLF(38:43)', nan, pre_HFLF(44:60)'];

dys_HRV = mean(pre_HRV_new(dysPD));
nondys_HRV = mean(pre_HRV_new(nondysPD));
healthy_HRV = mean(pre_HRV_new(healthy));

dys_HFLF = mean(pre_HFLF_new(dysPD));
nondys_HFLF = mean(pre_HFLF_new(nondysPD));
healthy_HFLF = mean(pre_HFLF_new(healthy));

% Getting brain region wise recording (Do this with only 19 channel values)

frontal = {'fp1','fp2','f3','f4','fz'};
temporo_occipital = {'f7','t3','t5','o1','o2','t6','t4','f8'};
centro_parietal = {'c3','cz','c4','p3','pz','p4'};

EEG_preproc = pop_loadset('C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\V1\1001\T1(11-25-07).set'); %Load just one EEG set file

% Getting the indices for the brain regions
for i=1:size(frontal,2)
    frontal_idx(i) = find(strcmpi({EEG_preproc.chanlocs.labels},frontal{i}));
end

for i=1:size(temporo_occipital,2)
    temporo_occipital_idx(i) = find(strcmpi({EEG_preproc.chanlocs.labels},temporo_occipital{i}));
end

for i=1:size(centro_parietal,2)
    centro_parietal_idx(i) = find(strcmpi({EEG_preproc.chanlocs.labels},centro_parietal{i}));
end


% PAC Features for subjects

count = 1;
for subidx = 1:length(dysPD)
    dys_normo_pre_F_T2(count) = mean(pac_normo{2,dysPD(subidx),1}(frontal_idx));
    dys_normo_F_T3(count) = mean(pac_normo{3,dysPD(subidx),3}(frontal_idx) - pac_normo{3,dysPD(subidx),1}(frontal_idx));
    dys_normo_TO_T3(count) = mean(pac_normo{3,dysPD(subidx),3}(temporo_occipital_idx) - pac_normo{3,dysPD(subidx),1}(temporo_occipital_idx));
    dys_tachy_pre_CP_T3(count) = mean(pac_tachy{3,dysPD(subidx),1}(centro_parietal_idx));
    count = count + 1;
end

dys_normo_pre_F_T2 = mean(dys_normo_pre_F_T2);
dys_normo_F_T3 = mean(dys_normo_F_T3);
dys_normo_TO_T3 = mean(dys_normo_TO_T3);
dys_tachy_pre_CP_T3 = mean(dys_tachy_pre_CP_T3);

count = 1;
for subidx = 1:length(nondysPD)
    nondys_normo_pre_F_T2(count) = mean(pac_normo{2,nondysPD(subidx),1}(frontal_idx));
    nondys_normo_F_T3(count) = mean(pac_normo{3,nondysPD(subidx),3}(frontal_idx) - pac_normo{3,nondysPD(subidx),1}(frontal_idx));
    nondys_normo_TO_T3(count) = mean(pac_normo{3,nondysPD(subidx),3}(temporo_occipital_idx) - pac_normo{3,nondysPD(subidx),1}(temporo_occipital_idx));
    nondys_tachy_pre_CP_T3(count) = mean(pac_tachy{3,nondysPD(subidx),1}(centro_parietal_idx));
    count = count + 1;
end

nondys_normo_pre_F_T2 = mean(nondys_normo_pre_F_T2);
nondys_normo_F_T3 = mean(nondys_normo_F_T3);
nondys_normo_TO_T3 = mean(nondys_normo_TO_T3);
nondys_tachy_pre_CP_T3 = mean(nondys_tachy_pre_CP_T3);

count = 1;
for subidx = 1:length(healthy)
    healthy_normo_pre_F_T2(count) = mean(pac_normo{2,healthy(subidx),1}(frontal_idx));
    healthy_normo_F_T3(count) = mean(pac_normo{3,healthy(subidx),3}(frontal_idx) - pac_normo{3,healthy(subidx),1}(frontal_idx));
    healthy_normo_TO_T3(count) = mean(pac_normo{3,healthy(subidx),3}(temporo_occipital_idx) - pac_normo{3,healthy(subidx),1}(temporo_occipital_idx));
    healthy_tachy_pre_CP_T3(count) = mean(pac_tachy{3,healthy(subidx),1}(centro_parietal_idx));
    count = count + 1;
end

healthy_normo_pre_F_T2 = mean(healthy_normo_pre_F_T2);
healthy_normo_F_T3 = mean(healthy_normo_F_T3);
healthy_normo_TO_T3 = mean(healthy_normo_TO_T3);
healthy_tachy_pre_CP_T3 = mean(healthy_tachy_pre_CP_T3);

% Spider plot

dys_pac = [dys_normo_pre_F_T2, dys_normo_F_T3, dys_normo_TO_T3, dys_tachy_pre_CP_T3];
nondys_pac = [nondys_normo_pre_F_T2, nondys_normo_F_T3, nondys_normo_TO_T3, nondys_tachy_pre_CP_T3];
healthy_pac = [healthy_normo_pre_F_T2, healthy_normo_F_T3, healthy_normo_TO_T3, healthy_tachy_pre_CP_T3];

dys_heart = [dys_HRV, dys_HFLF];
nondys_heart = [nondys_HRV, nondys_HFLF];
healthy_heart = [healthy_HRV, healthy_HFLF];

minlim_pac = min(min([dys_pac; nondys_pac; healthy_pac]))*1.1;
maxlim_pac = max(max([dys_pac; nondys_pac; healthy_pac]))*1.1;
minlim_heart = min(min([dys_heart; nondys_heart; healthy_heart]))*0.8;
maxlim_heart = max(max([dys_heart; nondys_heart; healthy_heart]))*1.5;

dys = [dys_pac, dys_heart];
nondys = [nondys_pac, nondys_heart];
healthy = [healthy_pac, healthy_heart];

s = spider_plot_class([dys; nondys; healthy]);

% Spider plot properties
s.AxesLimits = [minlim_pac minlim_pac minlim_pac minlim_pac 0 0;...
    maxlim_pac maxlim_pac maxlim_pac maxlim_pac 4 4]; % [min axes limits; max axes limits]
s.AxesLabels = {'PAC normo pre F T2','PAC normo F T3','PAC normo TO T3','PAC tachy pre CP T3', 'pre HRV', 'pre LF/HF ratio'};
s.FillOption = {'on','on','on'};
s.FillTransparency = [0.2,0.2,0.2];
s.LegendLabels = {'PD dyskinesia', 'PD non-dyskinesia', 'Healthy'};
s.AxesPrecision = 2;


%% Calculating PAC
%{
TO_labels = {'f7','t3','t5','o1','o2','t6','t4','f8'};

for i = 1:size(TO_labels,2)
     %Getting the eeg signal and bandpass filtering it
     eeg_sig = eeg_Healthy.data(find(strcmpi({eeg_Healthy.chanlocs.labels},TO_labels{i})),:);
     eeg_sig = bandpass(eeg_sig',[13,30],250);

     eegH = hilbert(eeg_sig)';
     eggH = hilbert(egg_sig_Healthy);

     eeg_amp = abs(eegH);
     egg_phase = angle(eggH);
     eeg_egg_coupled_signal = eeg_amp.*exp(1i*egg_phase);

     percentile_threshold = prctile(eeg_amp, 95);
     indices_top_5_percentile = find(eeg_amp >= percentile_threshold);

     %Top 5%ile PAC
     pac_numerator = abs(sum(eeg_egg_coupled_signal(indices_top_5_percentile)));
     pac_denominator = sqrt(sum(eeg_amp(indices_top_5_percentile).*eeg_amp(indices_top_5_percentile)))*sqrt(length(eeg_amp(indices_top_5_percentile)));
     %pac_top5percentile = pac_numerator/pac_denominator;
     pac_top5percentile = pac_numerator / sqrt(sum(eeg_amp(indices_top_5_percentile).*eeg_amp(indices_top_5_percentile)));

     %% Polar plot
     phase = angle(eeg_egg_coupled_signal(indices_top_5_percentile));
     magnitude = abs(eeg_egg_coupled_signal(indices_top_5_percentile));
     angle_PAC = angle(sum(eeg_egg_coupled_signal(indices_top_5_percentile)));
     denom = sqrt((eeg_amp(indices_top_5_percentile).*eeg_amp(indices_top_5_percentile)));
     magnitude = magnitude./denom;
     
     figure;
     polarplot([repelem(0,length(phase)); phase], ...
         [repelem(0,length(phase));magnitude],'b'); hold on;     
     polarplot([0;angle_PAC],[0;pac_top5percentile],'r','LineWidth',3);
     title(append('PD(',TO_labels{i},')'));
     %rlim([0 0.6]);
     %legend({'Vectors',repmat({''},length(phase)-1),'Mean Vector Length'})
end

%}