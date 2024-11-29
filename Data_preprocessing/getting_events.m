addpath(genpath('C:\Users\sanke\Downloads\Sanket\eeglab2023.1'))
visits = [1,3];
subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002 , 1029, 1050 not being considered
%subs_to_run = [1065];
set_path = 'C:\Users\sanke\Downloads\Sanket\preprocessed_EEG_files\';
mat_path = 'C:\Users\sanke\Downloads\Sanket\';


%EEG = pop_loadset('T1(13-28-42).set','C:\Users\sanke\Downloads\Sanket\preprocessed_EEG_files\V1\1065');
%eegplot(EEG.data)
events = [];

for visit = visits

    for subs = subs_to_run

        fileList = dir(append(set_path ,'V' , num2str(visit), '\' , num2str(subs) , '\*.set')); %getting the list of .set files in the directory
        
        %making a new directory for the power values to be saved
        %file_path = append(mat_path, 'EEG_power_files\V', num2str(visit), '\', num2str(subs));
        %mkdir(file_path);

        for files = 1:size(fileList,1)

            setfile = append(set_path ,'V' , num2str(visit), '\' , num2str(subs) , '\');

            %% Loading the dataset 

            EEG_preproc = pop_loadset(fileList(files,1).name,setfile);
            events(files,subs-1000,visit) = size(EEG_preproc.event,2);

        end
    end
end
