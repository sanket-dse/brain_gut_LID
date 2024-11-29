load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\hrv.mat')
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\lfhf.mat')
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\hf.mat')
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\lf.mat')

% 
% Getting post-pre features
post_hrv = hrv(:,:,2);
pre_hrv = hrv(:,:,1);
post_lfhf = lfhf(:,:,2);
pre_lfhf = lfhf(:,:,1);
%}
%{
%Getting pre features
hf = hf(:,:,1);
hrv = hrv(:,:,1);
lfhf = lfhf(:,:,1);
lf = lf(:,:,1);
%}
features = [pre_hrv post_hrv pre_lfhf post_lfhf];
% Turning them into a table
features_table = array2table(features,'VariableNames',{'pre_HRV_T1','pre_HRV_T2','pre_HRV_T3','pre_HRV_T4', ...
                                                       'post_HRV_T1','post_HRV_T2','post_HRV_T3','post_HRV_T4',...
                                                       'pre_HFLF_T1','pre_HFLF_T2','pre_HFLF_T3','pre_HFLF_T4', ...
                                                       'post_HFLF_T1','post_HFLF_T2','post_HFLF_T3','post_HFLF_T4'});

% Removing some subjects due to absence of all electrophysiological features
features_table([21,24,35,41],:) = [];

%% Converting the table to csv file
writetable(features_table,'heart_features.csv','Delimiter',',');