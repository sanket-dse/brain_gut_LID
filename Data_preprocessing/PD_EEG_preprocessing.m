%% Adding files to path and initializing some variables

addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1\'));
visits = [1,3];
subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002,1029 and 1050 not considered
%subs_to_run = [1007];
mat_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\mat_files\';
set_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\';
missing_channels = nan(4,67,3);

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

%% Preprocessing the EEG data 

for visit = visits

    for subs = subs_to_run

        fileList = dir(append(mat_path ,'V' , num2str(visit), '\' , num2str(subs) , '\*.mat')); %getting the list of .mat files in the directory
        
        %making a new directory for the preprocessed set files to be saved
        file_path = append(set_path, 'preprocessed_EEG_files\V', num2str(visit), '\', num2str(subs));
        mkdir(file_path);

        for files = 1:size(fileList,1)

            setfile = append(mat_path ,'V' , num2str(visit), '\' , num2str(subs) , '\' , fileList(files,1).name);

            %% Loading the .mat file into the EEG struct 

            EEG_raw = pop_importdata('dataformat','matlab','nbchan',0,'data',setfile,'srate',256,'pnts',0,'xmin',0);
             
            %% Creating the EEG struct
            EEG_raw.trials = 1;
            EEG_raw.nbchan = size(EEG_raw.data,1);
            EEG_raw.pnts = size(EEG_raw.data,2);
            EEG_raw.srate = 256; %Hz sampling rate
            EEG_raw.xmin = 0;
            EEG_raw.xmax = size(EEG_raw.data,2)/EEG_raw.srate;
            EEG_raw.times = linspace(EEG_raw.xmin,EEG_raw.xmax,EEG_raw.pnts);
            EEG_raw.etc = [];

            %reading channel locations 
            EEG_raw.chanlocs = chanlocs_new;
            
            %% Resampling the EEG_raw data to match that of EGG sampling rate

            EEG_raw = pop_resample(EEG_raw,250);

            %% Trimming some data at the end because it is just flat 

            end_rem_dur = 3; %seconds
            EEG_pop = pop_select(EEG_raw,"time",[EEG_raw.xmin EEG_raw.xmax-end_rem_dur]);
            fprintf('Data removal done');

            %% Filtering the data (FIR filter)

            EEG_filt = pop_eegfiltnew(EEG_pop, 'locutoff',  1, 'hicutoff',  45, 'filtorder', 9000, 'plotfreqz', 0); % this gives better results

            %% Remove bad channels using clean_channels function

            EEG_filt_rmchan = clean_channels(EEG_filt);
            fprintf('Bad channels removed');

            %% Performing ICA

            ch_list = 1:EEG_filt_rmchan.nbchan;

            EEG_filt_rmchan = pop_runica(EEG_filt_rmchan,'runica');
            EEG_filt_rmchan.icachansind = double(1:EEG_filt_rmchan.nbchan);
            EEG_filt_rmchan = iclabel(EEG_filt_rmchan);
    
            fprintf('ICA done\n')

            %% Remove bad components and reconstruct data
            threshold_signal = 0.05; % brain component with less than 5% confidence is removed
            cls = EEG_filt_rmchan.etc.ic_classification.ICLabel.classes;
            cls_score = EEG_filt_rmchan.etc.ic_classification.ICLabel.classifications;
            bad_comp = [];
            for cmp=1:size(EEG_filt_rmchan.icachansind,2)
                if cls_score(cmp,1)<threshold_signal
                    bad_comp = [bad_comp,cmp];
                end
            end
            EEG_ica = pop_subcomp(EEG_filt_rmchan, bad_comp, 0);
    
            fprintf('Bad components removed\n');

            %% Identifying bad portions of the data and reconstructing using ASR

            EEG_clean = clean_rawdata(EEG_ica,[],-1,-1,[],[],[]);
            fprintf('Bad portions of data identified and reconstructed \n');
            
            %% Calculating missing channels (temporary)
            task = str2num(fileList(files,1).name(2));
            missing_channels(task,subs-1000,visit) = 19 - size(EEG_clean.data,1);

            %% Interpolate bad channels
            EEG_interpol = pop_interp(EEG_clean, EEG_raw.chanlocs, 'spherical');
            fprintf('Bad channels interpolated\n');

            %% Re-reference the data
            EEG_reref = pop_reref(EEG_interpol,[]);
            fprintf('Data rereferenced\n');

            %% Save set
            %EEG_reref.setname = fileList(files,1).name(1:2);
            %EEG_reref.filepath = file_path;
            %EEG_reref.filename = append(fileList(files,1).name(1:end-4),'.set');
            %pop_saveset(EEG_reref,'savemode','resave');
    
            %fprintf('Saved dataset\n');
        end
    end
end

%% Calculating the missing channel percentage

missing_channels(:,:,2) = [];
missing_channels(:,[2,29,50],:) = [];
missing_channels(:,[21,24,35,41],:) = [];

missing_percentage = nanmean(nanmean(missing_channels,1),3); 


