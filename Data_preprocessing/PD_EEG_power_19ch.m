addpath(genpath('C:\Users\sanke\Downloads\Sanket_EEG\eeglab2023.1'));
addpath('C:\Users\sanke\Downloads\PD_EEG_EGG\fooof_mat-main\fooof_mat');
visits = [1,3];
subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002 , 1029, 1050 not being considered
%subs_to_run = [1038];
set_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\preprocessed_EEG_files\';
mat_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\';


%% Running the loop to get power

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
            EEG_data = EEG_preproc.data;

            %% Calculating the average power in beta band 

            fs = EEG_preproc.srate;
            nfft = 5*fs;
            window = hanning(nfft); %rectwin()
            noverlap = nfft/2;
            %FOOOF settings
            settings = struct();
            f_range = [1, 30];

            for c=1:size(EEG_data,1)
                [Pxx,w] =  pwelch(EEG_data(c,:),window,noverlap,nfft,fs);
                power_norm(c,:) = Pxx/sum(Pxx);     
            end
            
            %power_no_aperiodic_norm = power_no_aperiodic./ sum(power_no_aperiodic,2); %Normalizing the power spectral density

            task = str2num(fileList(files,1).name(2)); %get the task number from the filename
                      
            entire_ind = find(w>=0 & w<=30);
            b_ind = find(w>0.5 & w<4);
            %beta_max_pow(1,:) = max(power_no_aperiodic(:,b_ind),[],2);
            beta_tot_pow{task,subs-1000,visit} = sum(power_norm(:,b_ind),2);
            beta_mean_pow{task,subs-1000,visit} = mean(power_norm(:,b_ind),2);
            beta_relative_pow{task,subs-1000,visit} = sum(power_norm(:,b_ind),2)./sum(power_norm(:,entire_ind),2); %Normalizing w.r.t the entire 0-30Hz spectrum

        end
    end
end

save('delta_EEG_Power_19ch.mat','beta_relative_pow','beta_mean_pow','beta_tot_pow');

%
%% Plotting the EEG power 

%{
for visit = [1,3]

    for task = 1:4

        figure;
        
        count = 1;
        for sub = subs_to_run            
            subplot(8,8,count)
            if isempty(beta_relative_pow{task,sub-1000,visit})
                continue
            end

            topoplot(beta_relative_pow{task,sub-1000,visit},EEG_preproc.chanlocs);
            colorbar;
            caxis([min(beta_relative_pow{task,sub-1000,visit}),max(beta_relative_pow{task,sub-1000,visit})])
            title(['Subject ',num2str(sub-1000)])
            count = count + 1;
        end
        sgtitle(['For visit ',num2str(visit),' Task ',num2str(task)])

    end
end

%}
