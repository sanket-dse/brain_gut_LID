%% Loading some files and initializing variables

addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
raw_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\mat_files\';
preproc_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\';

missing_points = [];
visits = [1,3];
subs_to_run = 1001:1067;
% Removing the subject numbers who weren't considered
% PS - their numbers = 1002, 1022, 1025, 1029, 1037, 1043, 1050 
subs_to_run([2,22,25,29,37,43,50])=[];



%% Running the loop to get the points removed

for subs = subs_to_run
    
    for visit = visits
        
        %Getting the list of files in the directory
        fileList_set = dir(append(preproc_path ,'V' , num2str(visit), '\' , num2str(subs) , '\*.set'));
        fileList_mat = dir(append(raw_path ,'V' , num2str(visit), '\' , num2str(subs) , '\*.mat'));
        
        for files = 1:size(fileList_set,1)
            
            % Getting the task number
            task = str2num(fileList_mat(files,1).name(2));

            % load the raw EEG mat file
            load(append(raw_path,'\V',num2str(visit),'\',num2str(subs),'\',fileList_mat(files,1).name));
            
            if visit==1 
                if ~exist('data_V1')
                    raw_EEG = data;
                else
                    raw_EEG = data_V1;
                end

            else     
                if ~exist('data_V3')
                    raw_EEG = data;
                else
                    raw_EEG = data_V3;
                end
            end

            % load the preprocessed EEG set file
            EEG_preproc = pop_loadset(append(preproc_path,'\V',num2str(visit),'\',num2str(subs),'\',fileList_set(files,1).name)); %Load EEG set file
            preproc_EEG = EEG_preproc.data;
            
            % Getting the difference in points for the raw and preprocessed EEG
            diff = size(raw_EEG,2) - size(preproc_EEG,2) + 3*EEG_preproc.srate; %Compensating for the 3 seconds that were removed as flatline
            perc_diff = (diff/size(raw_EEG,2))*100;

            missing_points(task,subs-1000,visit) = perc_diff; 

            clear preproc_EEG
            clear -regexp ^data
            
        end
    end
end

%% Averaging the removed points

% Removing the subjects 
missing_points(:,[2,22,25,29,37,43,50],:) =[];
missing_points(:,:,2) = [];
missing_points(missing_points<0) = nan;

%Taking mean of the missing points
grand_mean = mean(mean(missing_points),3);
fprintf('Subjectwise mean of removed points : %f\n',mean(grand_mean));
fprintf('Subjectwise std of removed points : %f\n',std(grand_mean));