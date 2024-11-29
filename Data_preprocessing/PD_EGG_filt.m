addpath(genpath('/data3/Sanket/eeglab2023.1'));
egg_path = '/data4/NEURO_DrRathour/DATA/OPEN BCI/';
filt_path = '/data3/Sanket/';
visits = [3];
%subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002,1029 and 1050 not being run 
subs_to_run = [1047:1049,1051:1067]; %subs 1002,1029 and 1050 not being run 
for visit = visits
    
    for sub = subs_to_run

        if ~isfile(append(egg_path,'V' , num2str(visit),'/',num2str(sub),'.txt'))
            fprintf('%d.txt doesnt exist',sub)
            continue
        end
                   

        egg_raw_file = readtable(append(egg_path,'V' , num2str(visit),'/',num2str(sub),'.txt'));
        egg_sig = table2array(egg_raw_file(:,6:9));
        egg_sig = fillmissing(egg_sig, 'previous');
        fprintf('EGG file loaded\n')        
       
        % Filter signal
        fprintf('Being filtered\n')
        %s1_filt = bandpass(egg_sig,[0.03 0.07],250); %Normogastric range
        %s1_filt = bandpass(egg_sig,[0.0083 0.03],250); %Bradygastric range
        s1_filt = bandpass(egg_sig,[0.07 0.15],250); %Tachygastric range
        fprintf('FIltering done\n')
        
        %Create a directory
        mkdir(append(filt_path,'filtered_EGG_tachy_range/V',num2str(visit),'/'));
        save(append(filt_path,'filtered_EGG_tachy_range/V',num2str(visit),'/',num2str(sub),'.mat'),'s1_filt'); %saving the file
    end
end