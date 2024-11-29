% Adding the toolbox to the path
addpath(genpath('C:\Users\sanke\Downloads\PD_EEG_EGG\MediationToolbox-master'))
addpath(genpath('C:\Users\sanke\Downloads\PD_EEG_EGG\CanlabCore-master'))

% Importing the features 
features = readtable('C:\Users\sanke\Downloads\PD_models\Complete_Features_v3.csv');
heart_features = readtable('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\Demographics_heart_features.csv');

%% EGG power -> EEG power -> Dyskinesia Score (multilevel Mediation)
% Getting the variables
M = sqrt(features.eeg_power_pre_CP_T1(features.PD_Healthy==1));
M(:,2) = sqrt(features.eeg_power_pre_CP_T2(features.PD_Healthy==1));
M(:,3) = sqrt(features.eeg_power_pre_CP_T3(features.PD_Healthy==1));
M(:,4) = sqrt(features.eeg_power_pre_CP_T4(features.PD_Healthy==1));

Y = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,2) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,3) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,4) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);

X = features.tachy_power_pre_T1(features.PD_Healthy==1);
X(:,2) = features.tachy_power_pre_T2(features.PD_Healthy==1);
X(:,3) = features.tachy_power_pre_T3(features.PD_Healthy==1);
X(:,4) = features.tachy_power_pre_T4(features.PD_Healthy==1);

%Imputing values
%X = knnimpute(X);
%M = knnimpute(M);
%Y = knnimpute(Y);

Age = features.Age(features.PD_Healthy == 1)./80;
Gender_M = dummyvar(categorical(features.Gender(features.PD_Healthy == 1)));
Gender_M = Gender_M(:,2);

%%%%%% Run these three lines for single level, else skip
%X = reshape(X,1,144);
%Y = reshape(Y,1,144);
%M = reshape(M,1,144);

X_c = num2cell(X,1);
M_c = num2cell(M,1);
Y_c = num2cell(Y,1);

% Running the mediation analysis
[paths, stats] = mediation(X_c, Y_c, M_c,'verbose','boot','plots','names', {'Pre condition Tachy EGG power','Dyskinesia score ','Pre condition Centro-parietal EEG power'});

%% For PAC -> Heart feature -> Dyskinesia score
% Getting the variables
X = log(features.pac_normo_pre_TO_T1(features.PD_Healthy==1))./4;
X(:,2) = log(features.pac_normo_pre_TO_T2(features.PD_Healthy==1))./4;
X(:,3) = log(features.pac_normo_pre_TO_T3(features.PD_Healthy==1))./4;
X(:,4) = log(features.pac_normo_pre_TO_T4(features.PD_Healthy==1))./4;

Y = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,2) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,3) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,4) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);

M = log(heart_features.pre_HRV_T1)./3;
M(:,2) = log(heart_features.pre_HRV_T2)./3;
M(:,3) = log(heart_features.pre_HRV_T3)./3;
M(:,4) = log(heart_features.pre_HRV_T4)./3;
M(M==-inf) = nan;

X_c = num2cell(X,1);
M_c = num2cell(M,1);
Y_c = num2cell(Y,1);

% Running the mediation analysis
[paths, stats] = mediation(X_c, Y_c, M_c,'verbose','boot','plots','names', {'PAC brady pre CP','Dyskinesia score ','Pre HFLF'});
%% Multilevel mediation with custom code for random permutation test (EGG power -> EEG power -> Dyskinesia Score)

% Getting the variables
M = features.eeg_power_CP_T1(features.PD_Healthy==1);
M(:,2) = features.eeg_power_CP_T2(features.PD_Healthy==1);
M(:,3) = features.eeg_power_CP_T3(features.PD_Healthy==1);
M(:,4) = features.eeg_power_CP_T4(features.PD_Healthy==1);

Y = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,2) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,3) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,4) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);

X = features.tachy_power_T1(features.PD_Healthy==1);
X(:,2) = features.tachy_power_T2(features.PD_Healthy==1);
X(:,3) = features.tachy_power_T3(features.PD_Healthy==1);
X(:,4) = features.tachy_power_T4(features.PD_Healthy==1);

X_c = num2cell(X,1);
M_c = num2cell(M,1);
Y_c = num2cell(Y,1);

%Random permutation test
for nboots = 1:1000 
    Y_b = []; M_b = []; X_b = [];

    for i=1:3
        X_b(:,i) = X(randperm(size(X,1)),i);
        Y_b(:,i) = Y(randperm(size(Y,1)),i);
        M_b(:,i) = M(randperm(size(M,1)),i);
    end

    X_b = num2cell(X_b,1);
    M_b = num2cell(M_b,1);
    Y_b = num2cell(Y_b,1);

    [paths, stats] = mediation(X_b, Y_b, M_b);
    a(nboots) = stats.mean(1);
    b(nboots) = stats.mean(2);
    c_prime(nboots) = stats.mean(3);
    c(nboots) = stats.mean(4);
    ab(nboots) = stats.mean(5);

end

[paths, stats] = mediation(X_c, Y_c, M_c,'plots');
a_est = stats.mean(1);
b_est = stats.mean(2);
c_prime_est = stats.mean(3);
c_est = stats.mean(4);
ab_est = stats.mean(5);


%% EEG power -> PAC -> Dyskinesia Score (Multilevel Moderation)
% Although I am not sure if the methodology used is right?
% Note : We need to impute missing values to use moderation

% Getting the variables
X = features.eeg_power_CP_T1(features.PD_Healthy==1);
X(:,2) = features.eeg_power_CP_T2(features.PD_Healthy==1);
X(:,3) = features.eeg_power_CP_T3(features.PD_Healthy==1);
X(:,4) = features.eeg_power_CP_T4(features.PD_Healthy==1);

Y = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,2) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,3) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,4) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);

M = features.pac_tachy_F_T1(features.PD_Healthy==1);
M(:,2) = features.pac_tachy_F_T2(features.PD_Healthy==1);
M(:,3) = features.pac_tachy_F_T3(features.PD_Healthy==1);
M(:,4) = features.pac_tachy_F_T4(features.PD_Healthy==1);

%Imputing values
X = knnimpute(X);
M = knnimpute(M);
Y = knnimpute(Y);


X_c = num2cell(X,1);
M_c = num2cell(M,1);
Y_c = num2cell(Y,1);

% Running the mediation analysis
[paths, stats,statslow,statshigh] = moderation(X_c, Y_c, M_c,'verbose','plots','boot');

%% Single level moderation 

% Getting the variables
X = features.eeg_power_CP_T1(features.PD_Healthy==1);
X(:,2) = features.eeg_power_CP_T2(features.PD_Healthy==1);
X(:,3) = features.eeg_power_CP_T3(features.PD_Healthy==1);
X(:,4) = features.eeg_power_CP_T4(features.PD_Healthy==1);

Y = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,2) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,3) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);
Y(:,4) = sqrt(features.DyskinesiaScore(features.PD_Healthy==1)./8);

M = features.pac_normo_CP_T1(features.PD_Healthy==1);
M(:,2) = features.pac_normo_CP_T2(features.PD_Healthy==1);
M(:,3) = features.pac_normo_CP_T3(features.PD_Healthy==1);
M(:,4) = features.pac_normo_CP_T4(features.PD_Healthy==1);

%for i=1:4; X(:,i) = scale(X(:,i));end
%for i=1:4; M(:,i) = scale(M(:,i));end 

X = reshape(X,1,144);
Y = reshape(Y,1,144);
M = reshape(M,1,144);

X = (X - nanmean(X))./nanstd(X);
M = (M - nanmean(M))./nanstd(M);

% Running the bootstrap
b0 = [];
b1 = [];
b2 = [];
b3 = [];
for nboots = 1:1000 
    [X_b,idx] = datasample(X,length(X));
    Y_b = Y(idx);
    M_b = M(idx);
    lm_tbl = table(X_b',Y_b',M_b','VariableNames',{'X','Y','M'});
    int_model = fitlm(lm_tbl,'Y ~ X*M');
    b0(nboots) = int_model.Coefficients.Estimate(1);
    b1(nboots) = int_model.Coefficients.Estimate(2);
    b2(nboots) = int_model.Coefficients.Estimate(3);
    b3(nboots) = int_model.Coefficients.Estimate(4);
end
b0 = sort(b0);
b1 = sort(b1);
b2 = sort(b2);
b3 = sort(b3);

% Print the confidence intervals
fprintf('b0 : 2.5th percentile = %f       97.5th percentile = %f\n',b0(25),b0(975));
fprintf('b1 : 2.5th percentile = %f       97.5th percentile = %f\n',b1(25),b1(975));
fprintf('b2 : 2.5th percentile = %f       97.5th percentile = %f\n',b2(25),b2(975));
fprintf('b3 : 2.5th percentile = %f       97.5th percentile = %f\n',b3(25),b3(975));

%%Simple slopes procedure (For M = mean, mean+sd, mean-sd)
%For M = mean

M_level = nanmean(M);
x = min(X):0.1:max(X);
y = mean(b0) + (mean(b1) + mean(b3)*M_level).*x + M_level*mean(b2);
figure;plot(x,y);hold on;

M_level = nanmean(M) + nanstd(M);
y = mean(b0) + (mean(b1) + mean(b3)*M_level).*x + M_level*mean(b2);
plot(x,y);hold on;

M_level = nanmean(M) - nanstd(M);
y = mean(b0) + (mean(b1) + mean(b3)*M_level).*x + M_level*mean(b2);
plot(x,y);
xlabel('CP EEG power');
ylabel('Dyskinesia Score');
legend({'mean M','(mean + sd) M','(mean-sd) M'})