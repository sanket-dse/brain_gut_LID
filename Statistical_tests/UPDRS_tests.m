path1 = 'C:\Users\sanke\Downloads\PD_EEG_EGG\NURO Record.xlsx';
[~,sheets] = xlsfinfo(path1);
[~, ~, patient_record] = xlsread(path1,sheets{:,1});
[~, ~, UPDRS] = xlsread(path1,sheets{:,3});

%% Extracting all the data from patient records

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

%% Visualizing some stuff

x = ["Males" "Females"];
y = [length(find(subject_info(find(subject_info(:,2)==1),4)==77)) length(find(subject_info(find(subject_info(:,2)==1),4)==70))];
figure;
bar(x,y)
title('Gender(PD patients)');

y = [length(find(subject_info(find(subject_info(:,2)==0),4)==77)) length(find(subject_info(find(subject_info(:,2)==0),4)==70))];
figure;
bar(x,y)
title('Gender(healthy people)');
%% Extracting data from UPDRS score sheet

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
%all_scores(idx,2:end) = nan;

%Wilcoxon tests on all the scores
addpath(genpath('C:\Users\sanke\Downloads\PD_EEG_EGG\p-value_adjust'));
addpath('C:\Users\sanke\Downloads\PD_EEG_EGG\Violinplot-Matlab-master');

p_value = nan(size(all_scores,2)-1,1);

for i = 2:size(all_scores,2)
    if i==30
        p_value(i-1) = nan;
        continue
    end
    healthy = all_scores(find(subject_info(:,2)==0),i);
    patient = all_scores(find(subject_info(:,2)==1),i);
    p_value(i-1) = ranksum(healthy,patient);
end

adjusted_p = pval_adjust(p_value,'bonferroni');

arr = find(adjusted_p <= 0.05);
p_val = adjusted_p(arr);
sig_p = [arr,p_val];

% Creating a readable and referenceable sheet
UPDRS_new = UPDRS(:,[3:30,34:75]);

%% Visualizing differences in the scores

PD_Healthy = [repelem(1,34),repelem(0,11),1,1,0,0,0,1,0,0,0,0,0,1,repelem(0,7)];
healthy_idx = find(PD_Healthy == 0);
PD_idx = find(PD_Healthy == 1);

PD_scores = all_scores(PD_idx,71);
healthy_scores = all_scores(healthy_idx,71);
healthy_scores(end:length(PD_scores)) = nan;


figure;
boxplot([healthy_scores,PD_scores],'Labels',{'Healthy','PD'});

%% Calculating correlation between UPDRS scores and EEG EGG measures

load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\EEG_EGG_PAC_new.mat');
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\EEG_power.mat');

% Some pac values to remove 
rem_arr = [1,65,3; %task , sub id, visit
           2,53,3;
           1,43,3;
           3,66,1];
for i=1:length(rem_arr)
    pac_beta_top5{rem_arr(i,1),rem_arr(i,2),rem_arr(i,3)} = nan(1,3);
    pac_beta{rem_arr(i,1),rem_arr(i,2),rem_arr(i,3)} = nan(1,3);
end

% Extracting values from the pac structure for frontal,T1,V1
pac_frontal_T1_V1 = [];
for i=1:size(pac_beta_top5,2)
    if isempty(pac_beta{1,i,1})
        pac_frontal_T1_V1(i) = nan;
        continue
    end
    pac_frontal_T1_V1(i) = pac_beta{1,i,1}(1);
end
pac_frontal_T1_V1 = pac_frontal_T1_V1([1,3:28,30:49,51:67]); %Removing subs 1002,1029,1050
rho = corr(pac_frontal_T1_V1',all_scores(:,2),'rows','complete','type','Spearman');

% Calculating correlation with the UPDRS scores 

pac = [];

% getting PAC values
for i = 1:size(pac_beta_top5,2) % subjects
    for j = 1:size(pac_beta_top5,1) % tasks
        for k = [1,3] %visits
            for m=1:3 %regions
            
            if isempty(pac_beta{j,i,k})
                pac(j,i,k,m) = nan;
                continue
            end

            pac(j,i,k,m) = pac_beta{j,i,k}(m);
            end
        end
    end
end

pac = pac(:,[1,3:28,30:49,51:67],:,:); %Removing subs 1002,1029,1050

%Calculating correlations

addpath('C:\Users\sanke\Downloads\PD_EEG_EGG\findND'); %Adding the repo for find function for N dimensional matrices

rho_pac_V1 = [];
rho_pac_V3 = [];
pval_pac_V1 = [];
pval_pac_V3 = [];
for j = 1:size(pac_beta_top5,1) % tasks
    for m=1:3 %regions
        for i = 2:size(all_scores,2)
            [rho_pac_V1(i-1,j,m),pval_pac_V1(i-1,j,m)] = corr(pac(j,:,1,m)',all_scores(:,i),'rows','complete','type','Spearman');
            [rho_pac_V3(i-1,j,m),pval_pac_V3(i-1,j,m)]  = corr(pac(j,:,3,m)',all_scores(:,i),'rows','complete','type','Spearman');
        end
    end
end

%finding the locations of pvals < 0.05 (significant correlations) 
[i1,j1,k1] = findND(pval_pac_V1 < 0.05);
[i2,j2,k2] = findND(pval_pac_V3 < 0.05);

%Finding the value of significant correlations
for count=1:length(i1)
    corr_V1(count) = rho_pac_V1(i1(count),j1(count),k1(count));
    pval_V1(count) = pval_pac_V1(i1(count),j1(count),k1(count));
end

for count=1:length(i2)
    corr_V3(count) = rho_pac_V3(i2(count),j2(count),k2(count));
    pval_V3(count) = pval_pac_V3(i2(count),j2(count),k2(count));
end

%% Plotting the scores vs Pac values
%{
for count=1:length(i1)
    figure;
    scatter(pac(j1(count),:,1,k1(count))',all_scores(:,i1(count)));
end

for count=1:length(i2)
    figure;
    scatter(pac(j2(count),:,1,k2(count))',all_scores(:,i2(count)));
end
%}
