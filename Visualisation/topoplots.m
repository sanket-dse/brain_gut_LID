%% Load files and toolboxes

load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG\beta_PAC_normo_range_19ch.mat');
load('C:\Users\sanke\Downloads\PD_EEG_EGG\variables\beta_EEG\beta_EEG_Power_19ch.mat');
addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));

%% Getting chanlocs
% channel list 
labels = {'fp2','f4','c4','p4','fp1','f3','c3','p3','f8','t4','t6','o2','f7','t3','t5','o1','cz','pz','fz'};
%reading the chanlocs txt file 

chansinfo = readmatrix('C:\Users\sanke\Downloads\Sanket_EEG\chanlocs.txt');
temp = readtable('C:\Users\sanke\Downloads\Sanket_EEG\chanlocs.txt');
chanlabels = table2array(temp(:,1));
chansinfo = chansinfo(:,2:end);
chanlocs_new = [];
nchans = length(labels);

for i=1:nchans

    idx = find(strcmpi(chanlabels,labels{i}));
    chanlocs_new(:,i).X = chansinfo(idx,1);
    chanlocs_new(:,i).Y = chansinfo(idx,2);
    chanlocs_new(:,i).Z = chansinfo(idx,3);

    [sph_theta, sph_phi, chanlocs_new(:,i).sph_radius] = cart2sph(chanlocs_new(:,i).X,chanlocs_new(:,i).Y,chanlocs_new(:,i).Z);
    chanlocs_new(:,i).sph_phi = rad2deg(sph_phi);
    chanlocs_new(:,i).sph_theta = rad2deg(sph_theta);

    [chanlocs_new(:,i).urchanlocs] = i;
    [~,chanlocs_new(:,i).theta,chanlocs_new(:,i).radius] = sph2topo([chanlocs_new(:,i).urchanlocs,chanlocs_new(:,i).sph_phi,chanlocs_new(:,i).sph_theta]);
    chanlocs_new(:,i).labels = char(chanlabels(idx));

    chanlocs_new(:,i).theta = chanlocs_new(:,i).theta + 90 ; %this was done because the positions were rotated by 90 degrees


end
%% Averaging the values

rel_pow = [];
pac = [];
for subs = 1:67
    for visits = 1:3
        for tasks = 1:4

            if isempty(beta_relative_pow{tasks,subs,visits})
                rel_pow(:,subs,visits,tasks) = nan(19,1);
                continue
            end
            rel_pow(:,subs,visits,tasks) = beta_relative_pow{tasks,subs,visits};
    
            if isempty(pac_beta_top5{tasks,subs,visits})
                pac(:,subs,visits,tasks) = nan(19,1);
                continue
            end
            pac(:,subs,visits,tasks) = pac_beta_top5{tasks,subs,visits};
        end
    end
end

% Removing missing points
pac(:,:,2,:) = [];
pac(:,[2,29,50],:,:) = [];
rel_pow(:,:,2,:) = [];

%Taking the mean across tasks and visits and subjects
rel_pow2 = [];
pac2 = []
count = 1;
for tasks = 1:4
    for visits = 1:2
        rel_pow2(:,count) = nanmean(rel_pow(:,:,visits,tasks),2);
        pac2(:,count) = nanmean(pac(:,:,visits,tasks),2);
        count = count + 1;
    end
end

mean_pac = mean(pac2,2);
mean_rel_pow = mean(rel_pow2,2);

%% Topoplots

EEG.chanlocs = chanlocs_new;
% EEG power
figure;topoplot(mean_rel_pow,EEG.chanlocs,'maplimits',[min(mean_rel_pow),max(mean_rel_pow)]);
colorbar;
caxis([min(mean_rel_pow),max(mean_rel_pow)])
title('Beta EEG power')
% Broadband PAC
figure;topoplot(mean_pac,EEG.chanlocs,'maplimits',[min(mean_pac),max(mean_pac)]);
colorbar;
caxis([min(mean_pac),max(mean_pac)])
title('Broadband PAC')



    




