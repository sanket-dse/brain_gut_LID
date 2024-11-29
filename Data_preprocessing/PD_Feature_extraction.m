%% Loading some files and initializing variables
addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\EGG_power_v4.mat');
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG\beta_EEG_Power_19ch.mat');
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG_2nd_run\beta_PAC_normo_range_19ch_2.mat');
pac_normo = pac_beta_top5;
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG_2nd_run\beta_PAC_tachy_range_19ch_2.mat');
pac_tachy = pac_beta_top5;
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG_2nd_run\beta_PAC_brady_range_19ch_2.mat');
pac_brady = pac_beta_top5;
clear('pac_beta_top5','pac_beta');

subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002,1029 and 1050 not considered
PD_Healthy = [repelem(1,34),repelem(0,11),1,1,0,0,0,1,0,0,0,0,0,1,repelem(0,7)]; %1-> patient, 0-> healthy
visits = [1,3];
PD = find(PD_Healthy == 1 | PD_Healthy == 0); %Use this to toggle between only PD and all subjects

% Variables
brady_power_anova = [];
brady_power_pre_anova = [];
brady_power_post_anova = [];
pac_tachy_anova = [];
pac_tachy_pre_anova = [];
pac_tachy_post_anova = [];
pac_normo_anova = [];
pac_normo_pre_anova = [];
pac_normo_post_anova = [];
pac_brady_anova = [];
pac_brady_pre_anova = [];
pac_brady_post_anova = [];
normo_power_anova = [];
normo_power_pre_anova = [];
normo_power_post_anova = [];
eeg_power_anova = [];
eeg_power_pre_anova = [];
eeg_power_post_anova = [];
tachy_power_anova = [];
tachy_power_pre_anova = [];
tachy_power_post_anova = [];

% Some pac values to remove due to excessive zero padding
rem_arr = [1,65,3; %task , sub id, visit
           2,53,3;
           1,43,3;
           3,66,1];

for i=1:length(rem_arr)
    pac_normo{rem_arr(i,1),rem_arr(i,2),rem_arr(i,3)} = nan(1,19);
    pac_tachy{rem_arr(i,1),rem_arr(i,2),rem_arr(i,3)} = nan(1,19);
    pac_brady{rem_arr(i,1),rem_arr(i,2),rem_arr(i,3)} = nan(1,19);
end

% Converting all 0s to NaN since they were not calculated
brady_power(find(brady_power==0)) = nan;
normo_power(find(normo_power==0)) = nan;
tachy_power(find(tachy_power==0)) = nan;

%% Getting brain region wise recording (Do this with only 19 channel values)

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

% Getting the data brain regionwise

for visit = visits

    for subs = subs_to_run

        for task = 1:4

            if ~isempty(pac_brady{task,subs-1000,visit})
                pac_brady_new{task,subs-1000,visit}(1) = mean(pac_brady{task,subs-1000,visit}(frontal_idx));
                pac_brady_new{task,subs-1000,visit}(2) = mean(pac_brady{task,subs-1000,visit}(temporo_occipital_idx));
                pac_brady_new{task,subs-1000,visit}(3) = mean(pac_brady{task,subs-1000,visit}(centro_parietal_idx));
            end

            if ~isempty(beta_relative_pow{task,subs-1000,visit})
                beta_relative_pow_new{task,subs-1000,visit}(1) = mean(beta_relative_pow{task,subs-1000,visit}(frontal_idx));
                beta_relative_pow_new{task,subs-1000,visit}(2) = mean(beta_relative_pow{task,subs-1000,visit}(temporo_occipital_idx));
                beta_relative_pow_new{task,subs-1000,visit}(3) = mean(beta_relative_pow{task,subs-1000,visit}(centro_parietal_idx));
            end

            if ~isempty(pac_normo{task,subs-1000,visit})
                pac_normo_new{task,subs-1000,visit}(1) = mean(pac_normo{task,subs-1000,visit}(frontal_idx));
                pac_normo_new{task,subs-1000,visit}(2) = mean(pac_normo{task,subs-1000,visit}(temporo_occipital_idx));
                pac_normo_new{task,subs-1000,visit}(3) = mean(pac_normo{task,subs-1000,visit}(centro_parietal_idx));
            end

            if ~isempty(pac_tachy{task,subs-1000,visit})
                pac_tachy_new{task,subs-1000,visit}(1) = mean(pac_tachy{task,subs-1000,visit}(frontal_idx));
                pac_tachy_new{task,subs-1000,visit}(2) = mean(pac_tachy{task,subs-1000,visit}(temporo_occipital_idx));
                pac_tachy_new{task,subs-1000,visit}(3) = mean(pac_tachy{task,subs-1000,visit}(centro_parietal_idx));
            end      
        end
    end
end

% For ease of coding let's reassign
pac_brady = pac_brady_new;
pac_normo = pac_normo_new;
pac_tachy = pac_tachy_new;
beta_relative_pow = beta_relative_pow_new;

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

% Removing subs 1002 and 1050 (since they are not considered for analysis)

idx = find(subject_info(:,1)==1002 | subject_info(:,1)==1050);
subject_info(idx,:) = nan;
subject_info = rmmissing(subject_info);

% Creating a readable and referenceable sheet
UPDRS_new = UPDRS(:,[3:30,34:75]);

%% Loop to create the power values and PAC values array  

count=1;
for sub = subs_to_run(PD)
    
    %for visit = visits
        
        for task = 1:4
            
            for region = 1:3               
                
                %%%%%%%% For tachy pac ################
                if isempty(pac_tachy{task,sub-1000,1}) || isempty(pac_tachy{task,sub-1000,3}) 
                    pac_tachy_anova(count,region) = nan;
                else
                    pac_tachy_anova(count,region) = pac_tachy{task,sub-1000,3}(region) - pac_tachy{task,sub-1000,1}(region);  
                end  

                if isempty(pac_tachy{task,sub-1000,1})                    
                    pac_tachy_pre_anova(count,region) = nan;                 
                else                    
                    pac_tachy_pre_anova(count,region) = pac_tachy{task,sub-1000,1}(region);                    
                end  

                if isempty(pac_tachy{task,sub-1000,3})                 
                    pac_tachy_post_anova(count,region) = nan;
                else                     
                    pac_tachy_post_anova(count,region) = pac_tachy{task,sub-1000,3}(region); 
                end  
                
                %%%%%%%%% For normo pac %%%%%%%%%%%%%
                if isempty(pac_normo{task,sub-1000,3}) || isempty(pac_normo{task,sub-1000,1})
                    pac_normo_anova(count,region) = nan;                   
                else
                    pac_normo_anova(count,region) = pac_normo{task,sub-1000,3}(region) - pac_normo{task,sub-1000,1}(region);                   
                end  

                if isempty(pac_normo{task,sub-1000,1})
                    pac_normo_pre_anova(count,region) = nan;
                else                    
                    pac_normo_pre_anova(count,region) = pac_normo{task,sub-1000,1}(region);                   
                end  

                if isempty(pac_normo{task,sub-1000,3})
                    pac_normo_post_anova(count,region) = nan;                    
                else                    
                    pac_normo_post_anova(count,region) = pac_normo{task,sub-1000,3}(region);
                end  

                %%%%%%%%%%% For brady pac %%%%%%%%%%%
                if isempty(pac_brady{task,sub-1000,3}) || isempty(pac_brady{task,sub-1000,1})
                    pac_brady_anova(count,region) = nan;                    
                else
                    pac_brady_anova(count,region) = pac_brady{task,sub-1000,3}(region) - pac_brady{task,sub-1000,1}(region) ;                    
                end  

                if isempty(pac_brady{task,sub-1000,1})
                    pac_brady_pre_anova(count,region) = nan;                    
                else                    
                    pac_brady_pre_anova(count,region) = pac_brady{task,sub-1000,1}(region) ;                    
                end  

                if isempty(pac_brady{task,sub-1000,3})
                    pac_brady_post_anova(count,region) = nan;                    
                else                    
                    pac_brady_post_anova(count,region) = pac_brady{task,sub-1000,3}(region) ;
                end  

                %%%%%%% For EEG power %%%%%%%%%
                
                if isempty(beta_relative_pow{task,sub-1000,3}) || isempty(beta_relative_pow{task,sub-1000,1})
                    eeg_power_anova(count,region) = nan;                  
                else
                    eeg_power_anova(count,region) = beta_relative_pow{task,sub-1000,3}(region) - beta_relative_pow{task,sub-1000,1}(region);                    
                end

                if isempty(beta_relative_pow{task,sub-1000,1})
                    eeg_power_pre_anova(count,region) = nan;                    
                else
                    eeg_power_pre_anova(count,region) = beta_relative_pow{task,sub-1000,1}(region);                    
                end

                if isempty(beta_relative_pow{task,sub-1000,3})
                    eeg_power_post_anova(count,region) = nan;                   
                else
                    eeg_power_post_anova(count,region) = beta_relative_pow{task,sub-1000,3}(region) ;                    
                end
            end   
                         
            % Loop to get EGG power anova arrays
            brady_power_anova(count) = brady_power(task,sub-1000,3) - brady_power(task,sub-1000,1);
            brady_power_pre_anova(count) = brady_power(task,sub-1000,1);
            brady_power_post_anova(count) = brady_power(task,sub-1000,3);

            normo_power_anova(count) = normo_power(task,sub-1000,3) - normo_power(task,sub-1000,1);
            normo_power_pre_anova(count) = normo_power(task,sub-1000,1);
            normo_power_post_anova(count) = normo_power(task,sub-1000,3);

            tachy_power_anova(count) = tachy_power(task,sub-1000,3) - tachy_power(task,sub-1000,1);
            tachy_power_pre_anova(count) = tachy_power(task,sub-1000,1);
            tachy_power_post_anova(count) = tachy_power(task,sub-1000,3);
            
            count = count+1;
        end
end
%end
        
%% Getting the features  

features = [];

%%%%%%%% post and pre PAC features %%%%%%%%%%

% Brady range pac_top5
for region=1:3
    if length(PD) == 64
        features = vertcat(features,reshape(pac_brady_pre_anova(:,region),[4,64]));
        features = vertcat(features,reshape(pac_brady_post_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(pac_brady_anova(:,region),[8,38])); %For PD patients
    end
end

%Normogastric pac_top5
for region=1:3
    if length(PD) == 64
        features = vertcat(features,reshape(pac_normo_pre_anova(:,region),[4,64]));
         features = vertcat(features,reshape(pac_normo_post_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(pac_normo_anova(:,region),[8,38])); %For PD patients
    end
end

%Tachygastric pac_top5
for region=1:3
    if length(PD) == 64
        features = vertcat(features,reshape(pac_tachy_pre_anova(:,region),[4,64]));
        features = vertcat(features,reshape(pac_tachy_post_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(pac_tachy_post_anova(:,region),[8,38])); %For PD patients
    end
end

%%%%%%%% post-pre PAC features %%%%%%%%%%%
% Brady range pac_top5
for region=1:3
    if length(PD) == 64
        features = vertcat(features,reshape(pac_brady_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(pac_brady_anova(:,region),[8,38])); %For PD patients
    end
end

%Normogastric pac_top5
for region=1:3
    if length(PD) == 64
        features = vertcat(features,reshape(pac_normo_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(pac_normo_anova(:,region),[8,38])); %For PD patients
    end
end
%Tachygastric pac_top5
for region=1:3
    if length(PD) == 64
        features = vertcat(features,reshape(pac_tachy_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(pac_tachy_anova(:,region),[8,38])); %For PD patients
    end
end

%%%%%%%%%%%% post and pre EEG and EGG power %%%%%%%%%%

%Bradygastric EGG power
if length(PD)==64
    features = vertcat(features,reshape(brady_power_pre_anova,[4,64]));
    features = vertcat(features,reshape(brady_power_post_anova,[4,64]));
else
    features = vertcat(features,reshape(brady_power_anova,[8,38])); %For PD patients
end

%Normogastric EGG power
if length(PD)==64
    features = vertcat(features,reshape(normo_power_pre_anova,[4,64]));
    features = vertcat(features,reshape(normo_power_post_anova,[4,64]));
else
    features = vertcat(features,reshape(normo_power_anova,[8,38])); %For PD patients
end

% Tachygastric EGG power
if length(PD)==64
    features = vertcat(features,reshape(tachy_power_pre_anova,[4,64]));
    features = vertcat(features,reshape(tachy_power_post_anova,[4,64]));
else
    features = vertcat(features,reshape(tachy_power_anova,[8,38])); %For PD patients
end

%Beta EEG Power
for region=1:3
    if length(PD)==64
        features = vertcat(features,reshape(eeg_power_pre_anova(:,region),[4,64]));
        features = vertcat(features,reshape(eeg_power_post_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(eeg_power_anova(:,region),[8,38]));
    end
end

%%%%% post-pre EEG and EGG power %%%%%%%%%%%

%Bradygastric EGG power
if length(PD)==64
    features = vertcat(features,reshape(brady_power_anova,[4,64]));    
else
    features = vertcat(features,reshape(brady_power_anova,[8,38])); %For PD patients
end

%Normogastric EGG power
if length(PD)==64
    features = vertcat(features,reshape(normo_power_anova,[4,64]));
else
    features = vertcat(features,reshape(normo_power_anova,[8,38])); %For PD patients
end

% Tachygastric EGG power
if length(PD)==64
    features = vertcat(features,reshape(tachy_power_anova,[4,64]));    
else
    features = vertcat(features,reshape(tachy_power_anova,[8,38])); %For PD patients
end

%Beta EEG Power
for region=1:3
    if length(PD)==64
        features = vertcat(features,reshape(eeg_power_anova(:,region),[4,64]));
    else
        features = vertcat(features,reshape(eeg_power_anova(:,region),[8,38]));
    end
end

% Creating a table of features



features_table = array2table(features','VariableNames',{'pac_brady_pre_F_T1','pac_brady_pre_F_T2','pac_brady_pre_F_T3','pac_brady_pre_F_T4', ...
                                                        'pac_brady_post_F_T1','pac_brady_post_F_T2','pac_brady_post_F_T3','pac_brady_post_F_T4', ...
                                                        'pac_brady_pre_TO_T1','pac_brady_pre_TO_T2','pac_brady_pre_TO_T3','pac_brady_pre_TO_T4', ...
                                                        'pac_brady_post_TO_T1','pac_brady_post_TO_T2','pac_brady_post_TO_T3','pac_brady_post_TO_T4', ...
                                                        'pac_brady_pre_CP_T1','pac_brady_pre_CP_T2','pac_brady_pre_CP_T3','pac_brady_pre_CP_T4', ...
                                                        'pac_brady_post_CP_T1','pac_brady_post_CP_T2','pac_brady_post_CP_T3','pac_brady_post_CP_T4', ...
                                                        'pac_normo_pre_F_T1','pac_normo_pre_F_T2','pac_normo_pre_F_T3','pac_normo_pre_F_T4', ...
                                                        'pac_normo_post_F_T1','pac_normo_post_F_T2','pac_normo_post_F_T3','pac_normo_post_F_T4', ...
                                                        'pac_normo_pre_TO_T1','pac_normo_pre_TO_T2','pac_normo_pre_TO_T3','pac_normo_pre_TO_T4', ...
                                                        'pac_normo_post_TO_T1','pac_normo_post_TO_T2','pac_normo_post_TO_T3','pac_normo_post_TO_T4', ...
                                                        'pac_normo_pre_CP_T1','pac_normo_pre_CP_T2','pac_normo_pre_CP_T3','pac_normo_pre_CP_T4', ...
                                                        'pac_normo_post_CP_T1','pac_normo_post_CP_T2','pac_normo_post_CP_T3','pac_normo_post_CP_T4', ...
                                                        'pac_tachy_pre_F_T1','pac_tachy_pre_F_T2','pac_tachy_pre_F_T3','pac_tachy_pre_F_T4', ...
                                                        'pac_tachy_post_F_T1','pac_tachy_post_F_T2','pac_tachy_post_F_T3','pac_tachy_post_F_T4', ...
                                                        'pac_tachy_pre_TO_T1','pac_tachy_pre_TO_T2','pac_tachy_pre_TO_T3','pac_tachy_pre_TO_T4', ...
                                                        'pac_tachy_post_TO_T1','pac_tachy_post_TO_T2','pac_tachy_post_TO_T3','pac_tachy_post_TO_T4', ...
                                                        'pac_tachy_pre_CP_T1','pac_tachy_pre_CP_T2','pac_tachy_pre_CP_T3','pac_tachy_pre_CP_T4', ...
                                                        'pac_tachy_post_CP_T1','pac_tachy_post_CP_T2','pac_tachy_post_CP_T3','pac_tachy_post_CP_T4', ...
                                                        ...
                                                        'pac_brady_F_T1','pac_brady_F_T2','pac_brady_F_T3','pac_brady_F_T4', ...                                                        
                                                        'pac_brady_TO_T1','pac_brady_TO_T2','pac_brady_TO_T3','pac_brady_TO_T4', ...                                                        
                                                        'pac_brady_CP_T1','pac_brady_CP_T2','pac_brady_CP_T3','pac_brady_CP_T4', ...                                                        
                                                        'pac_normo_F_T1','pac_normo_F_T2','pac_normo_F_T3','pac_normo_F_T4', ...                                                        
                                                        'pac_normo_TO_T1','pac_normo_TO_T2','pac_normo_TO_T3','pac_normo_TO_T4', ...                                                        
                                                        'pac_normo_CP_T1','pac_normo_CP_T2','pac_normo_CP_T3','pac_normo_CP_T4', ...                                                        
                                                        'pac_tachy_F_T1','pac_tachy_F_T2','pac_tachy_F_T3','pac_tachy_F_T4', ...                                                        
                                                        'pac_tachy_TO_T1','pac_tachy_TO_T2','pac_tachy_TO_T3','pac_tachy_TO_T4', ...                                                        
                                                        'pac_tachy_CP_T1','pac_tachy_CP_T2','pac_tachy_CP_T3','pac_tachy_CP_T4', ...
                                                        ...
                                                        'brady_power_pre_T1','brady_power_pre_T2','brady_power_pre_T3','brady_power_pre_T4', ...
                                                        'brady_power_post_T1','brady_power_post_T2','brady_power_post_T3','brady_power_post_T4', ...
                                                        'normo_power_pre_T1','normo_power_pre_T2','normo_power_pre_T3','normo_power_pre_T4', ...
                                                        'normo_power_post_T1','normo_power_post_T2','normo_power_post_T3','normo_power_post_T4', ...
                                                        'tachy_power_pre_T1','tachy_power_pre_T2','tachy_power_pre_T3','tachy_power_pre_T4', ...
                                                        'tachy_power_post_T1','tachy_power_post_T2','tachy_power_post_T3','tachy_power_post_T4', ...
                                                        'eeg_power_pre_F_T1','eeg_power_pre_F_T2','eeg_power_pre_F_T3','eeg_power_pre_F_T4', ...
                                                        'eeg_power_post_F_T1','eeg_power_post_F_T2','eeg_power_post_F_T3','eeg_power_post_F_T4', ...
                                                        'eeg_power_pre_TO_T1','eeg_power_pre_TO_T2','eeg_power_pre_TO_T3','eeg_power_pre_TO_T4', ...
                                                        'eeg_power_post_TO_T1','eeg_power_post_TO_T2','eeg_power_post_TO_T3','eeg_power_post_TO_T4', ...
                                                        'eeg_power_pre_CP_T1','eeg_power_pre_CP_T2','eeg_power_pre_CP_T3','eeg_power_pre_CP_T4', ...
                                                        'eeg_power_post_CP_T1','eeg_power_post_CP_T2','eeg_power_post_CP_T3','eeg_power_post_CP_T4',...
                                                        ...
                                                        'brady_power_T1','brady_power_T2','brady_power_T3','brady_power_T4', ...                                                        
                                                        'normo_power_T1','normo_power_T2','normo_power_T3','normo_power_T4', ...                                                        
                                                        'tachy_power_T1','tachy_power_T2','tachy_power_T3','tachy_power_T4', ...                                                        
                                                        'eeg_power_F_T1','eeg_power_F_T2','eeg_power_F_T3','eeg_power_F_T4', ...                                                        
                                                        'eeg_power_TO_T1','eeg_power_TO_T2','eeg_power_TO_T3','eeg_power_TO_T4', ...                                                        
                                                        'eeg_power_CP_T1','eeg_power_CP_T2','eeg_power_CP_T3','eeg_power_CP_T4'});




% For gender
Gender = {};
for i=1:size(subject_info,1)
    if subject_info(i,4) == 77
        Gender{i} = 'Male';
    else
        Gender{i} = 'Female';
    end
end

%Adding age and gender as properties
Gender = categorical(Gender);
features_table.Age = subject_info(:,3);
features_table.Gender = Gender';
% Adding PD Healthy remark to the table
features_table.PD_Healthy = PD_Healthy';

%Adding the Part total scores
part4 = all_scores(PD,71);
features_table.PartIVScore = part4;
part3 = all_scores(:,64);
features_table.PartIIIScore = part3;
part2 = all_scores(:,29);
features_table.PartIIScore = part2;
part1 = all_scores(:,15);
features_table.PartIScore = part1;
total_score = part1+part2+part3+part4;
features_table.TotalScore = total_score;
score = all_scores(PD,54)+all_scores(PD,55); %Part 3.15
features_table.Part3_15 = score;
score = all_scores(PD,56)+all_scores(PD,57); %Part 3.16
features_table.Part3_16 = score;
score = all_scores(PD,58)+all_scores(PD,59)+all_scores(PD,60)+all_scores(PD,61)+all_scores(PD,62); %Part 3.17
features_table.Part3_17 = score;
score = all_scores(PD,63); %Part 3.18
features_table.Part3_18 = score;

% Summing up dyskinesia related scores 4.1 & 4.2
part4 = all_scores(PD,65) + all_scores(PD,66);
part4(isnan(part4)) = 0;
features_table.DyskinesiaScore = part4;

% Removing some subjects due to absence of all electrophysiological features
features_table([21,24,35,41],:) = [];
%features_table([21,24],:) = []; %For PD patients

%% Converting the table to csv file
writetable(features_table,'C:/Users/sanke/Downloads/Complete_Features_v3.csv','Delimiter',',');
