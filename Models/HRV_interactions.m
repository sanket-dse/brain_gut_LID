%% Loading initial files
features = readtable('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\Selected_features_v3_linear.csv');
heart_measures = readtable('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\Demographics_heart_features.csv');
heart_measures = heart_measures(:,[2:end-1]); %Trimming the index and PD healthy columns
%colnames = features.Properties.VariableNames;

%% Demographics test 
demo_heart = readtable('C:\Users\sanke\Downloads\PD_models\pre_heart_features.csv');
demo_age = readtable('C:\Users\sanke\Downloads\PD_models\Complete_Features_v3.csv');
demo_heart = [demo_heart , array2table(demo_age.PD_Healthy)];

% For age 
PD = demo_age.Age(demo_age.PD_Healthy==1);
Healthy = demo_age.Age(demo_age.PD_Healthy==0);
fprintf('Mean age of PD (in years) : %f \n', mean(PD));
fprintf('Std age of PD (in years) : %f \n', std(PD));
fprintf('Mean age of Healthy (in years) : %f \n', mean(Healthy));
fprintf('Std age of Healthy (in years) : %f \n', std(Healthy));
[p h stats] = ranksum(PD,Healthy)

% For Dyskinesia score 
PD = demo_age.DyskinesiaScore(demo_age.PD_Healthy==1);
Healthy = demo_age.DyskinesiaScore(demo_age.PD_Healthy==0);
fprintf('Mean score of PD : %f \n', median(PD));
fprintf('Std score of PD : %f \n', std(PD));
fprintf('Mean score of Healthy : %f \n', median(Healthy));
fprintf('Std score of Healthy : %f \n', std(Healthy));
[p h stats] = ranksum(PD,Healthy)

% For heart measures
HRV = mean(demo_heart(:,[1,2,3,4]),2);
HRV.Properties.VariableNames = {'HRV'};
HFLF = mean(demo_heart(:,[13,14,15,16]),2);
HFLF.Properties.VariableNames = {'HFLF'};
heart_table = [HRV, HFLF, table(demo_heart.Var1,'VariableNames',{'PD_Healthy'})];

fprintf('Mean HFLF ratio of PD  : %f \n', nanmean(heart_table.HFLF(heart_table.PD_Healthy==1)));
fprintf('Std HFLF ratio of PD  : %f \n', nanstd(heart_table.HFLF(heart_table.PD_Healthy==1)));
fprintf('Mean HFLF ratio of Healthy : %f \n', nanmean(heart_table.HFLF(heart_table.PD_Healthy==0)));
fprintf('Std HFLF ratio of Healthy : %f \n', nanstd(heart_table.HFLF(heart_table.PD_Healthy==0)));
[p h stats] = ranksum(heart_table.HFLF(heart_table.PD_Healthy==1),heart_table.HFLF(heart_table.PD_Healthy==0))

fprintf('Mean HRV of PD  : %f \n', nanmean(heart_table.HRV(heart_table.PD_Healthy==1)));
fprintf('Std HRV of PD  : %f \n', nanstd(heart_table.HRV(heart_table.PD_Healthy==1)));
fprintf('Mean HRV of Healthy : %f \n', nanmean(heart_table.HRV(heart_table.PD_Healthy==0)));
fprintf('Std HRV of Healthy : %f \n', nanstd(heart_table.HRV(heart_table.PD_Healthy==0)));
[p h stats] = ranksum(heart_table.HRV(heart_table.PD_Healthy==1),heart_table.HRV(heart_table.PD_Healthy==0))

%% Creating a robust linear model

%Base model
lm_model = fitlm(features,'RobustOpts','on')

%% Creating models adding heart measures and demographics as covariates in the model 

%Creating the demographics and heart measures table
pre_HRV = table2array(mean(heart_measures(:,[1,2,3,4]),2)); %Taking mean across the tasks
pre_HFLF = table2array(mean(heart_measures(:,[5,6,7,8]),2));
demo_table = table(pre_HRV,pre_HFLF);
demo_table.Age = heart_measures.Age;
demo_table.Gender_M = heart_measures.Gender_M;
lm_table = [features, demo_table];
colnames = {'pac_normo_F_T3', 'pac_normo_pre_F_T2'}; %Only significant features

%Only demographics
modelSpec = '';
for i = 1:size(colnames,2)
    modelSpec = append(modelSpec, '+', colnames{i});
end
modelSpec = modelSpec(2:end);
modelSpec = append('y~',modelSpec,'+Age+Gender_M');

demo_lm_model = fitlm(lm_table,modelSpec,'RobustOpts','on')

% Demographics and heart measures table
modelSpec = '';
for i = 1:size(colnames,2)
    modelSpec = append(modelSpec, '+', colnames{i});
end
modelSpec = modelSpec(2:end);
modelSpec = append('y~',modelSpec,'+Age+Gender_M+pre_HRV + pre_HFLF');

demo_heart_lm_model = fitlm(lm_table,modelSpec,'RobustOpts','on')

%% Adding interaction of heart features with the base model

meas = 'pre_HFLF'; %Change the heart measure from here
heart_meas = lm_table(:,[meas]);
%colnames = features.Properties.VariableNames;


%Creating the model specification

modelSpec = '';
for i = 1:size(colnames,2)
    modelSpec = append(modelSpec, '+', colnames{i});
end

modelSpec = append(modelSpec, '+', meas);
for i=1:size(colnames,2)
    modelSpec = append(modelSpec, '+', colnames{i}, ':', meas);
end

modelSpec = modelSpec(2:end);
modelSpec = append('y~',modelSpec);

%Creating the interaction model table and the robust linear model
new_lm_table = [features , heart_meas];
new_lm_model = fitlm(new_lm_table,modelSpec,'RobustOpts','on')
anova_mdl = anova(new_lm_model);

%% Check if any interactions are significant and plot the interactions
figure;plotInteraction(new_lm_model,'pre_HRV','pac_brady_pre_TO_T2','effects')
title('Interaction of pre HRV with pac brady TO T1');
%xlabel('pac brady pre CP T4');
%ylabel('Dyskinesia score');
