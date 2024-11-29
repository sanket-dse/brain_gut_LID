%% Loading some files 
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

%% Initializing variables
subs_to_run = subject_info(:,1)';
PD_Healthy = subject_info(:,2)'; %1-> patient, 0-> healthy

Dyskinesia = nan(1,length(subject_info)) ;
for i = 1:length(subject_info)
    if subject_info(i,5) > 0
        Dyskinesia(i) = 1;
    else
        Dyskinesia(i) = 0 ;
    end
end

visits = [1,3];
PD = find(PD_Healthy == 1 | PD_Healthy == 0); %Use this to toggle between only PD and all subjects

% Variables
brady_power_anova = [];
pac_tachy_anova = [];
pac_normo_anova = [];
pac_brady_anova = [];
normo_power_anova = [];
eeg_power_anova = [];
tachy_power_anova = [];

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

%% Loop to create the power values and PAC values array  
% by subtracting the post med values by pre med values
count=1;
for sub = subs_to_run(PD)
    
    for visit = visits
        
        for task = 1:4
            
            for region = 1:3
                
                %if isempty(beta_tot_pow{task,sub-1000,visit}) || beta_relative_pow{task,sub-1000,visit}(region) < 0
                
               
                if isempty(pac_tachy{task,sub-1000,visit}) 
                    pac_tachy_anova(count,region) = nan;
                else
                    pac_tachy_anova(count,region) = pac_tachy{task,sub-1000,visit}(region);                    
                end  

                if isempty(pac_normo{task,sub-1000,visit})
                    pac_normo_anova(count,region) = nan;
                else
                    pac_normo_anova(count,region) = pac_normo{task,sub-1000,visit}(region);                    
                end  

                if isempty(pac_brady{task,sub-1000,visit})
                    pac_brady_anova(count,region) = nan;
                else
                    pac_brady_anova(count,region) = pac_brady{task,sub-1000,visit}(region) ;                   
                end  

                if isempty(beta_relative_pow{task,sub-1000,visit})
                    eeg_power_anova(count,region) = nan;
                else
                    eeg_power_anova(count,region) = beta_relative_pow{task,sub-1000,visit}(region);
                end
            end      
                         
      
            % Loop to get EGG power anova arrays
            brady_power_anova(count) = brady_power(task,sub-1000,visit);
            normo_power_anova(count) = normo_power(task,sub-1000,visit);
            tachy_power_anova(count) = tachy_power(task,sub-1000,visit);
            
            count = count+1;
        end
    end
end
        

%% Creating the factor arrays for ANOVA

% For PD/Healthy
PD_Healthy_anova = [];

for sub = PD_Healthy
    PD_Healthy_anova = [PD_Healthy_anova,repelem(sub,8)];
end

% For Dyskinesia
Dyskinesia_anova = [];

for sub = Dyskinesia
    Dyskinesia_anova = [Dyskinesia_anova,repelem(sub,8)];
end

% For scores 

scores_anova = [];

for subs=1:size(all_scores,1)
    scores_anova = [scores_anova,repelem(all_scores(subs,64),8)]; %Part-III total score
end

% For visit
visit_anova = [repmat([repelem(1,4),repelem(3,4)],1,60)];

%For task
task_anova = [repmat([repmat([1,2,3,4],1,2)],1,60)];


%% Fitting rmanova model for PAC (All regions and EGG ranges combined)

%Creating the between factors and repeated measures table
meas_reshape = [];
for region = 1:3
    meas = log(pac_normo_anova(:,region)); %Change the variable from here
    meas_reshape_normo = reshape(meas,[8,60])';
    meas = log(pac_brady_anova(:,region));
    meas_reshape_brady = reshape(meas,[8,60])';
    meas = log(pac_tachy_anova(:,region));
    meas_reshape_tachy = reshape(meas,[8,60])';
    meas_reshape = [meas_reshape, meas_reshape_normo, meas_reshape_brady, meas_reshape_tachy];
end

% Creating the Dyskinesia and non Dyskinesia factors
PD_Healthy_cohort = {};
for i=1:length(PD_Healthy)
    if PD_Healthy(i)==1 
        PD_Healthy_cohort{i} = 'PD';
    else
        PD_Healthy_cohort{i} = 'Healthy';
    end
end

Dyskinesia_cohort = {};
for i=1:length(Dyskinesia)
    if Dyskinesia(i)==1 
        Dyskinesia_cohort{i} = 'Dys';
    else
        Dyskinesia_cohort{i} = 'Nondys';
    end
end

%{
Subject = [repelem({'Patient'},34),repelem({'Healthy'},11),{'Patient','Patient','Healthy','Healthy','Healthy','Patient','Healthy','Healthy', ...
    'Healthy','Healthy','Healthy','Patient'},repelem({'Healthy'},7)];

PD_Healthy = [repelem({'Patient'},34),repelem({'Healthy'},11),{'Patient','Patient','Healthy','Healthy','Healthy','Patient','Healthy','Healthy', ...
    'Healthy','Healthy','Healthy','Patient'},repelem({'Healthy'},7)];

Dyskinesia = [repelem({'Dys'},9),'Nondys',repelem({'Dys'},10),'Nondys',repelem({'Dys'},7),{'Nondys','Dys','Dys',...
    'Nondys','Dys','Nondys'},repelem({'Nondys'},11),{'Dys','Dys','Nondys','Nondys','Nondys','Dys'},repelem({'Nondys'},5),'Dys',...
    repelem({'Nondys'},7)];
Dyskinesia{51} = 'Nondys';
Dyskinesia{63} = 'Dys';
%}


Age = subject_info(:,3);
Gender = {};
for i=1:size(subject_info,1)
    if subject_info(i,4) == 77
        Gender{i} = 'Male';
    else
        Gender{i} = 'Female';
    end
end


between = table(PD_Healthy_cohort',Dyskinesia_cohort',Gender',Age, ...
    meas_reshape(:,1),meas_reshape(:,2),meas_reshape(:,3),meas_reshape(:,4),meas_reshape(:,5), meas_reshape(:,6),meas_reshape(:,7),meas_reshape(:,8), ...
    meas_reshape(:,9),meas_reshape(:,10),meas_reshape(:,11),meas_reshape(:,12),meas_reshape(:,13), meas_reshape(:,14),meas_reshape(:,15),meas_reshape(:,16), ...
    meas_reshape(:,17),meas_reshape(:,18),meas_reshape(:,19),meas_reshape(:,20),meas_reshape(:,21), meas_reshape(:,22),meas_reshape(:,23),meas_reshape(:,24), ...
    meas_reshape(:,25),meas_reshape(:,26),meas_reshape(:,27),meas_reshape(:,28),meas_reshape(:,29), meas_reshape(:,30),meas_reshape(:,31),meas_reshape(:,32), ...
    meas_reshape(:,33),meas_reshape(:,34),meas_reshape(:,35),meas_reshape(:,36),meas_reshape(:,37), meas_reshape(:,38),meas_reshape(:,39),meas_reshape(:,40), ...
    meas_reshape(:,41),meas_reshape(:,42),meas_reshape(:,43),meas_reshape(:,44),meas_reshape(:,45), meas_reshape(:,46),meas_reshape(:,47),meas_reshape(:,48), ...
    meas_reshape(:,49),meas_reshape(:,50),meas_reshape(:,51),meas_reshape(:,52),meas_reshape(:,53), meas_reshape(:,54),meas_reshape(:,55),meas_reshape(:,56), ...
    meas_reshape(:,57),meas_reshape(:,58),meas_reshape(:,59),meas_reshape(:,60),meas_reshape(:,61), meas_reshape(:,62),meas_reshape(:,63),meas_reshape(:,64), ...
    meas_reshape(:,65),meas_reshape(:,66),meas_reshape(:,67),meas_reshape(:,68),meas_reshape(:,69), meas_reshape(:,70),meas_reshape(:,71),meas_reshape(:,72), ...
    'VariableNames', {'PD_Healthy','Dyskinesia','Gender','Age', ...
    'F_Pre_T1_normo','F_Pre_T2_normo','F_Pre_T3_normo','F_Pre_T4_normo','F_Post_T1_normo','F_Post_T2_normo','F_Post_T3_normo','F_Post_T4_normo', ...
    'F_Pre_T1_brady','F_Pre_T2_brady','F_Pre_T3_brady','F_Pre_T4_brady','F_Post_T1_brady','F_Post_T2_brady','F_Post_T3_brady','F_Post_T4_brady', ...
    'F_Pre_T1_tachy','F_Pre_T2_tachy','F_Pre_T3_tachy','F_Pre_T4_tachy','F_Post_T1_tachy','F_Post_T2_tachy','F_Post_T3_tachy','F_Post_T4_tachy', ...
    'TO_Pre_T1_normo','TO_Pre_T2_normo','TO_Pre_T3_normo','TO_Pre_T4_normo','TO_Post_T1_normo','TO_Post_T2_normo','TO_Post_T3_normo','TO_Post_T4_normo', ...
    'TO_Pre_T1_brady','TO_Pre_T2_brady','TO_Pre_T3_brady','TO_Pre_T4_brady','TO_Post_T1_brady','TO_Post_T2_brady','TO_Post_T3_brady','TO_Post_T4_brady', ...
    'TO_Pre_T1_tachy','TO_Pre_T2_tachy','TO_Pre_T3_tachy','TO_Pre_T4_tachy','TO_Post_T1_tachy','TO_Post_T2_tachy','TO_Post_T3_tachy','TO_Post_T4_tachy', ...   
    'CP_Pre_T1_normo','CP_Pre_T2_normo','CP_Pre_T3_normo','CP_Pre_T4_normo','CP_Post_T1_normo','CP_Post_T2_normo','CP_Post_T3_normo','CP_Post_T4_normo', ...
    'CP_Pre_T1_brady','CP_Pre_T2_brady','CP_Pre_T3_brady','CP_Pre_T4_brady','CP_Post_T1_brady','CP_Post_T2_brady','CP_Post_T3_brady','CP_Post_T4_brady', ...
    'CP_Pre_T1_tachy','CP_Pre_T2_tachy','CP_Pre_T3_tachy','CP_Pre_T4_tachy','CP_Post_T1_tachy','CP_Post_T2_tachy','CP_Post_T3_tachy','CP_Post_T4_tachy'});

%Creating the within factors table
visit = [repmat([repelem({'Pre'},4),repelem({'Post'},4)],1,3*3)]';
task = [repmat([{'T1'},{'T2'},{'T3'},{'T4'}],1,2*3*3)]';
egg_range = [repmat([repelem({'Normogastric'},8),repelem({'Bradygastric'},8),repelem({'Tachygastric'},8)],1,3)]';
brain_region = [repelem({'Frontal'},24), repelem({'Temporo-occipital'},24), repelem({'Centro-parietal'},24)]';
within = table(visit,task,egg_range,brain_region,'VariableNames',{'Pre_post_med','Task','EGG_range','Brain_region'});
%within.Pre_post_med = categorical(within.Pre_post_med);
%within.Task = categorical(within.Task);
%within.Pre_post_med_Task = within.Pre_post_med.*within.Task;

%Fitting the model
rm = fitrm(between,'F_Pre_T1_normo-CP_Post_T4_tachy ~ Dyskinesia','WithinDesign',within,'WithinModel', ...
    "Pre_post_med + Task + EGG_range + Brain_region + Task:Pre_post_med + Brain_region:Pre_post_med + EGG_range:Pre_post_med + Brain_region:Task + EGG_range:Task + Brain_region:Task:Pre_post_med + EGG_range:Task:Pre_post_med");
ranovatbl = ranova(rm,'WithinModel', ...
   'Pre_post_med + Task + EGG_range + Brain_region + Task:Pre_post_med + Brain_region:Pre_post_med + EGG_range:Pre_post_med + Brain_region:Task + EGG_range:Task + Brain_region:Task:Pre_post_med + EGG_range:Task:Pre_post_med')

%rm = fitrm(between,'F_Pre_T1_normo-CP_Post_T4_tachy ~ Dyskinesia ','WithinDesign',within,'WithinModel', ...
%    "Pre_post_med + Task + EGG_range + Brain_region + Task:Pre_post_med + Brain_region:Pre_post_med + EGG_range:Pre_post_med + Brain_region:Task + EGG_range:Task + Brain_region:Task:Pre_post_med + EGG_range:Task:Pre_post_med");
%ranovatbl = ranova(rm,'WithinModel', ...
%    'Pre_post_med + Task + EGG_range + Brain_region + Task:Pre_post_med + Brain_region:Pre_post_med + EGG_range:Pre_post_med + Brain_region:Task + EGG_range:Task + Brain_region:Task:Pre_post_med + EGG_range:Task:Pre_post_med')

%Posthoc test
Mrm1 = multcompare(rm,'Pre_post_med','By','Task')
Mrm3 = multcompare(rm,'PD_Healthy','By','Task')
Mrm2 = multcompare(rm,'Pre_post_med','By','Dyskinesia')
Mrm4 = multcompare(rm,'Dyskinesia')
Mrm5 = multcompare(rm,'Pre_post_med')
Mrm6 = multcompare(rm,'Task')
Mrm7 = multcompare(rm,'EGG_range')
Mrm8 = multcompare(rm,'Pre_post_med','By','EGG_range')
Mrm9 = multcompare(rm,'PD_Healthy','By','EGG_range')
Mrm10 = multcompare(rm,'PD_Healthy','By','Brain_region')
Mrm11 = multcompare(rm,'Task','By','Brain_region')
Mrm12 = multcompare(rm,'Brain_region')

%% Plotting the main effects of task 

pac = [reshape(pac_brady_anova,480*3,1); reshape(pac_normo_anova,480*3,1); reshape(pac_tachy_anova,480*3,1)];
cohort = repmat(PD_Healthy_anova',9,1);
pre_post_med = repmat(visit_anova',9,1);
cog_task = repmat(task_anova',9,1);
brain_region = repmat([repelem({'F'},480), repelem({'TO'},480), repelem({'CP'},480)]',3,1);
egg_range = [repelem({'brady'},480*3), repelem({'normo'},480*3), repelem({'tachy'},480*3)]';
plot_table = table(pac,cohort,pre_post_med,cog_task,brain_region,egg_range,'VariableNames',{'PAC','PD_Healthy','Visit','Task','Brain_region','EGG_range'});

% Extracting the required values and plotting them
meas = plot_table([strcmpi(plot_table.EGG_range,'normo') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==1 & plot_table.Visit==1],:);
mean_val = nanmean(reshape(meas.PAC,4,36)');
std_val = nanstd(reshape(meas.PAC,4,36)');
figure;errorbar([1,2,3,4],mean_val,std_val,'o');
figure; boxplot(reshape(meas.PAC,4,36)','Labels',{'1','2','3','4'});
title('Frontal normogastric PAC for PD and pre-med');

meas = plot_table([strcmpi(plot_table.EGG_range,'normo') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==1 & plot_table.Visit==3],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,36)','Labels',{'1','2','3','4'});
title('Frontal normogastric PAC for PD and post-med');

meas = plot_table([strcmpi(plot_table.EGG_range,'normo') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==0 & plot_table.Visit==1],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,24)','Labels',{'1','2','3','4'});
title('Frontal normogastric PAC for healthy and pre-med');

meas = plot_table([strcmpi(plot_table.EGG_range,'normo') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==0 & plot_table.Visit==3],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,24)','Labels',{'1','2','3','4'});
title('Frontal normogastric PAC for healthy and post-med');

meas = plot_table([strcmpi(plot_table.EGG_range,'tachy') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==1 & plot_table.Visit==1],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,36)','Labels',{'1','2','3','4'});

meas = plot_table([strcmpi(plot_table.EGG_range,'tachy') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==1 & plot_table.Visit==3],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,36)','Labels',{'1','2','3','4'});

meas = plot_table([strcmpi(plot_table.EGG_range,'tachy') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==0 & plot_table.Visit==1],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,24)','Labels',{'1','2','3','4'});

meas = plot_table([strcmpi(plot_table.EGG_range,'tachy') & strcmpi(plot_table.Brain_region,'F') & plot_table.PD_Healthy==0 & plot_table.Visit==3],:);
figure; maineffectsplot(meas.PAC,meas.Task)
figure; boxplot(reshape(meas.PAC,4,24)','Labels',{'1','2','3','4'});

%% Fitting RMANOVA model with EEG power

%Creating the between factors and repeated measures table
meas_reshape = [];
for region = 1:3
    meas = sqrt(eeg_power_anova(:,region)); %Change the variable from here
    meas_reshape = [meas_reshape, reshape(meas,[8,60])'];
end

Age = subject_info(:,3);
Gender = {};
for i=1:size(subject_info,1)
    if subject_info(i,4) == 77
        Gender{i} = 'Male';
    else
        Gender{i} = 'Female';
    end
end


between = table(PD_Healthy_cohort',Dyskinesia_cohort',Gender',Age, ...
    meas_reshape(:,1),meas_reshape(:,2),meas_reshape(:,3),meas_reshape(:,4),meas_reshape(:,5), meas_reshape(:,6),meas_reshape(:,7),meas_reshape(:,8), ...
    meas_reshape(:,9),meas_reshape(:,10),meas_reshape(:,11),meas_reshape(:,12),meas_reshape(:,13), meas_reshape(:,14),meas_reshape(:,15),meas_reshape(:,16), ...
    meas_reshape(:,17),meas_reshape(:,18),meas_reshape(:,19),meas_reshape(:,20),meas_reshape(:,21), meas_reshape(:,22),meas_reshape(:,23),meas_reshape(:,24), ...
    'VariableNames', {'PD_Healthy','Dyskinesia','Gender','Age', ...
    'F_Pre_T1','F_Pre_T2','F_Pre_T3','F_Pre_T4','F_Post_T1','F_Post_T2','F_Post_T3','F_Post_T4', ...
    'TO_Pre_T1','TO_Pre_T2','TO_Pre_T3','TO_Pre_T4_','TO_Post_T1','TO_Post_T2','TO_Post_T3','TO_Post_T4', ...
    'CP_Pre_T1','CP_Pre_T2','CP_Pre_T3','CP_Pre_T4','CP_Post_T1','CP_Post_T2','CP_Post_T3','CP_Post_T4'});

%Creating the within factors table
visit = [repmat([repelem({'Pre'},4),repelem({'Post'},4)],1,3)]';
task = [repmat([{'T1'},{'T2'},{'T3'},{'T4'}],1,2*3)]';
%egg_range = [repmat([repelem({'Normogastric'},8),repelem({'Bradygastric'},8),repelem({'Tachygastric'},8)],1,3)]';
brain_region = [repelem({'Frontal'},8), repelem({'Temporo-occipital'},8), repelem({'Centro-parietal'},8)]';
within = table(visit,task,brain_region,'VariableNames',{'Pre_post_med','Task','Brain_region'});
%within.Pre_post_med = categorical(within.Pre_post_med);
%within.Task = categorical(within.Task);
%within.Pre_post_med_Task = within.Pre_post_med.*within.Task;

%Fitting the model
rm = fitrm(between,'F_Pre_T1-CP_Post_T4 ~ Dyskinesia','WithinDesign',within,'WithinModel', ...
    "Pre_post_med + Task + Brain_region + Task:Pre_post_med + Brain_region:Pre_post_med + Brain_region:Task + Brain_region:Task:Pre_post_med");
ranovatbl = ranova(rm,'WithinModel', ...
    'Pre_post_med + Task + Brain_region + Task:Pre_post_med + Brain_region:Pre_post_med + Brain_region:Task + Brain_region:Task:Pre_post_med')

%Posthoc test
Mrm1 = multcompare(rm,'Pre_post_med','By','Task')
Mrm3 = multcompare(rm,'Task','By','PD_Healthy')
Mrm2 = multcompare(rm,'Pre_post_med','By','PD_Healthy')
Mrm4 = multcompare(rm,'Dyskinesia')
Mrm5 = multcompare(rm,'Pre_post_med')
Mrm6 = multcompare(rm,'Task')
Mrm10 = multcompare(rm,'PD_Healthy','By','Brain_region')
Mrm11 = multcompare(rm,'Task','By','Brain_region')
Mrm11 = multcompare(rm,'Brain_region','By','Task')
Mrm12 = multcompare(rm,'Brain_region')
Mrm13 =  multcompare(rm,'Task','By','Brain_region')

%{
%For 3-way interactions (PD_Healthy:Task:Brain_region)
eeg_power = reshape(eeg_power_anova,1,512*3);
cohort = repmat(PD_Healthy_anova,1,3);
cog_task = repmat(task_anova,1,3);
eeg_region = [repelem({'Frontal'},1,512), repelem({'Temporo-occipital'},1,512), repelem({'Centro-parietal'},1,512)];
plot_table = table(eeg_power',cohort',cog_task',eeg_region','VariableNames',{'EEG_power','PD_Healthy','Task','Brain_region'});

%For frontal region 
idx = find(strcmp(plot_table.Brain_region,'Frontal'));
figure;interactionplot(plot_table.EEG_power(idx),[plot_table.PD_Healthy(idx), plot_table.Task(idx)],"varnames",{'PD Healthy','Task'});
sgtitle('For Frontal region in EEG power')

%For temporo-occipital region 
idx = find(strcmp(plot_table.Brain_region,'Temporo-occipital'));
figure;interactionplot(plot_table.EEG_power(idx),[plot_table.PD_Healthy(idx), plot_table.Task(idx)],"varnames",{'PD Healthy','Task'});
sgtitle('For Temporo-occipital region in EEG power')

%For centro-parietal region 
idx = find(strcmp(plot_table.Brain_region,'Centro-parietal'));
figure;interactionplot(plot_table.EEG_power(idx),[plot_table.PD_Healthy(idx), plot_table.Task(idx)],"varnames",{'PD Healthy','Task'});
sgtitle('For Centro-parietal region in EEG power')
%}

%% Fitting RMANOVA model with EGG power

%Creating the between factors and repeated measures table
meas_reshape = [];
meas = normo_power_anova; %Change the variable from here
meas_reshape_normo = reshape(meas,[8,60])';
meas = brady_power_anova; %Change the variable from here
meas_reshape_brady = reshape(meas,[8,60])';
meas = tachy_power_anova; %Change the variable from here
meas_reshape_tachy = reshape(meas,[8,60])';
meas_reshape = [meas_reshape_normo, meas_reshape_brady, meas_reshape_tachy];

Age = subject_info(:,3);
Gender = {};
for i=1:size(subject_info,1)
    if subject_info(i,4) == 77
        Gender{i} = 'Male';
    else
        Gender{i} = 'Female';
    end
end


between = table(PD_Healthy_cohort',Dyskinesia_cohort',Gender',Age, ...
    meas_reshape(:,1),meas_reshape(:,2),meas_reshape(:,3),meas_reshape(:,4),meas_reshape(:,5), meas_reshape(:,6),meas_reshape(:,7),meas_reshape(:,8), ...
    meas_reshape(:,9),meas_reshape(:,10),meas_reshape(:,11),meas_reshape(:,12),meas_reshape(:,13), meas_reshape(:,14),meas_reshape(:,15),meas_reshape(:,16), ...
    meas_reshape(:,17),meas_reshape(:,18),meas_reshape(:,19),meas_reshape(:,20),meas_reshape(:,21), meas_reshape(:,22),meas_reshape(:,23),meas_reshape(:,24), ...
    'VariableNames', {'PD_Healthy','Dyskinesia','Gender','Age', ...
    'normo_Pre_T1','normo_Pre_T2','normo_Pre_T3','normo_Pre_T4','normo_Post_T1','normo_Post_T2','normo_Post_T3','normo_Post_T4', ...
    'brady_Pre_T1','brady_Pre_T2','brady_Pre_T3','brady_Pre_T4_','brady_Post_T1','brady_Post_T2','brady_Post_T3','brady_Post_T4', ...
    'tachy_Pre_T1','tachy_Pre_T2','tachy_Pre_T3','tachy_Pre_T4','tachy_Post_T1','tachy_Post_T2','tachy_Post_T3','tachy_Post_T4'});

%Creating the within factors table
visit = [repmat([repelem({'Pre'},4),repelem({'Post'},4)],1,3)]';
task = [repmat([{'T1'},{'T2'},{'T3'},{'T4'}],1,2*3)]';
egg_range = [repelem({'Normogastric'},8),repelem({'Bradygastric'},8),repelem({'Tachygastric'},8)]';
%brain_region = [repelem({'Frontal'},8), repelem({'Temporo-occipital'},8), repelem({'Centro-parietal'},8)]';
within = table(visit,task,egg_range,'VariableNames',{'Pre_post_med','Task','EGG_range'});
%within.Pre_post_med = categorical(within.Pre_post_med);
%within.Task = categorical(within.Task);
%within.Pre_post_med_Task = within.Pre_post_med.*within.Task;

%Fitting the model
rm = fitrm(between,'normo_Pre_T1-tachy_Post_T4 ~ Dyskinesia','WithinDesign',within,'WithinModel', ...
    "Pre_post_med + Task + EGG_range + Task:Pre_post_med + EGG_range:Pre_post_med + EGG_range:Task + EGG_range:Task:Pre_post_med");
ranovatbl = ranova(rm,'WithinModel', ...
    'Pre_post_med + Task + EGG_range + Task:Pre_post_med + EGG_range:Pre_post_med + EGG_range:Task + EGG_range:Task:Pre_post_med')

%Posthoc test
Mrm1 = multcompare(rm,'Pre_post_med','By','Task')
Mrm3 = multcompare(rm,'PD_Healthy','By','Task')
Mrm2 = multcompare(rm,'Pre_post_med','By','PD_Healthy')
Mrm4 = multcompare(rm,'PD_Healthy')
Mrm5 = multcompare(rm,'Pre_post_med')
Mrm6 = multcompare(rm,'Task')
Mrm10 = multcompare(rm,'PD_Healthy','By','EGG_range')
Mrm11 = multcompare(rm,'Task','By','EGG_range')
Mrm12 = multcompare(rm,'EGG_range','By','Task')
Mrm10 = multcompare(rm,'Pre_post_med','By','EGG_range')

%{
%For 3-way interactions (PD_Healthy:Pre_post_med:EGG_range)
egg_power = [normo_power_anova, brady_power_anova , tachy_power_anova];
cohort = repmat(PD_Healthy_anova,1,3);
med = repmat(visit_anova,1,3);
range = [repelem({'Normogastric'},1,512), repelem({'Bradygastric'},1,512), repelem({'Tachygastric'},1,512)];
plot_table = table(egg_power',cohort',med',range','VariableNames',{'EGG_power','PD_Healthy','Pre_post_med','EGG_ranges'});

%For normogastric range 
idx = find(strcmp(plot_table.EGG_ranges,'Normogastric'));
figure;interactionplot(plot_table.EGG_power(idx),[plot_table.PD_Healthy(idx), plot_table.Pre_post_med(idx)],"varnames",{'PD Healthy','Pre post med'});
sgtitle('For normogastric range in EGG power')

%For bradygastric range 
idx = find(strcmp(plot_table.EGG_ranges,'Bradygastric'));
figure;interactionplot(plot_table.EGG_power(idx),[plot_table.PD_Healthy(idx), plot_table.Pre_post_med(idx)],"varnames",{'PD Healthy','Pre post med'});
sgtitle('For bradygastric range in EGG power')

%For tachygastric range 
idx = find(strcmp(plot_table.EGG_ranges,'Tachygastric'));
figure;interactionplot(plot_table.EGG_power(idx),[plot_table.PD_Healthy(idx), plot_table.Pre_post_med(idx)],"varnames",{'PD Healthy','Pre post med'});
sgtitle('For tachygastric range in EGG power')
%}

%% Independence of PAC and EEG, EGG features

% Spearman correlations between EEG power & EGG power with PAC values

%For normogastric PAC
pval = [];
fprintf('Normogastric PAC')
pac = reshape(pac_normo_anova,1,512*3);
eeg = reshape(eeg_power_anova,1,512*3);
egg = repmat(normo_power_anova,1,3);
[rho, pval(1)] = corr(pac',eeg', 'type','Spearman','Rows','complete')
[rho, pval(2)] = corr(pac',egg', 'type','Spearman','Rows','complete')

%For tachygastric PAC
fprintf('Tachygastric PAC')
pac = reshape(pac_tachy_anova,1,512*3);
eeg = reshape(eeg_power_anova,1,512*3);
egg = repmat(tachy_power_anova,1,3);
[rho, pval(3)] = corr(pac',eeg', 'type','Spearman','Rows','complete')
[rho, pval(4)] = corr(pac',egg', 'type','Spearman','Rows','complete')

%For bradygastric PAC
fprintf('Bradygastric PAC')
pac = reshape(pac_brady_anova,1,512*3);
eeg = reshape(eeg_power_anova,1,512*3);
egg = repmat(brady_power_anova,1,3);
[rho, pval(5)] = corr(pac',eeg', 'type','Spearman','Rows','complete')
[rho, pval(6)] = corr(pac',egg', 'type','Spearman','Rows','complete')

%Multiple comparison correction
[fdr,q,priori,R2] = mafdr(pval,'Method','polynomial','Showplot',true);

%% Fitting rmanova model for PAC (Reserve)
%{
%Creating the between factors and repeated measures table
meas_reshape = [];
region = 3;
meas = log(pac_normo_anova(:,region)); %Change the variable from here
meas_reshape_normo = reshape(meas,[8,60])';
meas = log(pac_brady_anova(:,region));
meas_reshape_brady = reshape(meas,[8,60])';
meas = log(pac_tachy_anova(:,region));
meas_reshape_tachy = reshape(meas,[8,60])';
meas_reshape = [meas_reshape, meas_reshape_normo, meas_reshape_brady, meas_reshape_tachy];


%Subject = [repelem({'Patient'},34),repelem({'Healthy'},11),{'Patient','Patient','Healthy','Healthy','Healthy','Patient','Healthy','Healthy', ...
%    'Healthy','Healthy','Healthy','Patient'},repelem({'Healthy'},7)];

%Subject = [repelem({'Patient'},9),'Healthy',repelem({'Patient'},10),'Healthy',repelem({'Patient'},7),{'Healthy','Patient','Patient',...
%    'Healthy','Patient','Healthy'},repelem({'Healthy'},11),{'Patient','Patient','Healthy','Healthy','Healthy','Patient'},repelem({'Healthy'},5),'Patient',...
%    repelem({'Healthy'},7)];
%Subject{51} = 'Healthy';
%Subject{63} = 'Patient';



Age = subject_info(:,3);
Gender = {};
for i=1:size(subject_info,1)
    if subject_info(i,4) == 77
        Gender{i} = 'Male';
    else
        Gender{i} = 'Female';
    end
end


between = table(Subject',Gender',Age, ...
    meas_reshape(:,1),meas_reshape(:,2),meas_reshape(:,3),meas_reshape(:,4),meas_reshape(:,5), meas_reshape(:,6),meas_reshape(:,7),meas_reshape(:,8), ...
    meas_reshape(:,9),meas_reshape(:,10),meas_reshape(:,11),meas_reshape(:,12),meas_reshape(:,13), meas_reshape(:,14),meas_reshape(:,15),meas_reshape(:,16), ...
    meas_reshape(:,17),meas_reshape(:,18),meas_reshape(:,19),meas_reshape(:,20),meas_reshape(:,21), meas_reshape(:,22),meas_reshape(:,23),meas_reshape(:,24), ...
    'VariableNames', {'PD_Healthy','Gender','Age', ...
    'Pre_T1_normo','Pre_T2_normo','Pre_T3_normo','Pre_T4_normo','Post_T1_normo','Post_T2_normo','Post_T3_normo','Post_T4_normo', ...
    'Pre_T1_brady','Pre_T2_brady','Pre_T3_brady','Pre_T4_brady','Post_T1_brady','Post_T2_brady','Post_T3_brady','Post_T4_brady', ...
    'Pre_T1_tachy','Pre_T2_tachy','Pre_T3_tachy','Pre_T4_tachy','Post_T1_tachy','Post_T2_tachy','Post_T3_tachy','Post_T4_tachy'});

%Creating the within factors table
visit = [repmat([repelem({'Pre'},4),repelem({'Post'},4)],1,3)]';
task = [repmat([{'T1'},{'T2'},{'T3'},{'T4'}],1,2*3)]';
egg_range = [repelem({'Normogastric'},8),repelem({'Bradygastric'},8),repelem({'Tachygastric'},8)]';
%brain_region = [repelem({'Frontal'},24), repelem({'Temporo-occipital'},24), repelem({'Centro-parietal'},24)]';
within = table(visit,task,egg_range,'VariableNames',{'Pre_post_med','Task','EGG_range'});

%Fitting the model
rm = fitrm(between,'Pre_T1_normo-Post_T4_tachy ~ PD_Healthy','WithinDesign',within,'WithinModel', ...
    "Pre_post_med + Task + EGG_range + Task:Pre_post_med + EGG_range:Pre_post_med + EGG_range:Task + EGG_range:Task:Pre_post_med");
ranovatbl = ranova(rm,'WithinModel', ...
    'Pre_post_med + Task + EGG_range + Task:Pre_post_med + EGG_range:Pre_post_med + EGG_range:Task + EGG_range:Task:Pre_post_med')

%Posthoc test
Mrm1 = multcompare(rm,'Pre_post_med','By','Task')
Mrm3 = multcompare(rm,'PD_Healthy','By','Task')
Mrm2 = multcompare(rm,'Pre_post_med','By','PD_Healthy')
Mrm4 = multcompare(rm,'PD_Healthy')
Mrm5 = multcompare(rm,'Pre_post_med')
Mrm6 = multcompare(rm,'Task')
Mrm7 = multcompare(rm,'EGG_range')
Mrm8 = multcompare(rm,'Pre_post_med','By','EGG_range')
Mrm9 = multcompare(rm,'PD_Healthy','By','EGG_range')
%}