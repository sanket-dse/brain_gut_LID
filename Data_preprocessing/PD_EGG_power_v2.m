%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PS : The output cell has the rows as tasks and the column is the
% subject. Inside each cell, the first element is frontal, 2nd element is
% temporo-occipital and 3rd element is centro-parietal. The third dimension
% is the visit (1 and 3).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tachy_path = '/data3/Sanket/EGG_tachy_taskwise/';
normo_path = '/data3/Sanket/EGG_normo_taskwise/';
brady_path = '/data3/Sanket/EGG_brady_taskwise/';
visits = [1,3];
%subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002 , 1029, 1050 not being considered
subs_to_run = [1001:1067];
%subs_to_run = [1001];

tachy_power = [];
normo_power = [];
brady_power = [];
%% Running the loop to get power

for visit = visits

    for subs = subs_to_run

        for tasks = 1:4
            
            %PSD estimation parameters
            fs_egg = 250;
            nfft = 60*fs_egg;
            window = hanning(nfft); %rectwin()
            noverlap = nfft/2;
            
            %Check whether file exists or not for normogastric
            if ~isfile(append(normo_path,'V' , num2str(visit),'/',num2str(subs), ...
                    '/Task',num2str(tasks),'.mat'))
                fprintf('%d .mat doesnt exist\n',subs)
                normo_power(tasks,subs-1000,visit) = nan;
            else
                %Loading the egg signal
                load(append(normo_path,'V' , num2str(visit),'/',num2str(subs), ...
                    '/Task',num2str(tasks),'.mat'));
                %Zero padding the signal if the length is less than the
                %window size
                if length(egg_section) < 60*fs_egg
                    egg_section(length(egg_section):60*fs_egg) = 0; 
                end
                %PSD estimation 
                [Pxx1,w1] =  pwelch(egg_section,window,noverlap,nfft,fs_egg);
                power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

                range_idx = find(w1>0.03 & w1<0.07); %Normogastric range
                normo_power(tasks,subs-1000,visit) = sum(power_norm(range_idx,:))./sum(power_norm); %Normogastric range

            end

            %Check whether file exists or not for tachygastric
            if ~isfile(append(tachy_path,'V' , num2str(visit),'/',num2str(subs), ...
                    '/Task',num2str(tasks),'.mat'))
                fprintf('%d .mat doesnt exist\n',subs)
                tachy_power(tasks,subs-1000,visit) = nan;
            else
                %Loading the egg signal
                load(append(tachy_path,'V' , num2str(visit),'/',num2str(subs), ...
                    '/Task',num2str(tasks),'.mat'));
                %Zero padding the signal if the length is less than the
                %window size
                if length(egg_section) < 60*fs_egg
                    egg_section(length(egg_section):60*fs_egg) = 0; 
                end
                %PSD estimation 
                [Pxx1,w1] =  pwelch(egg_section,window,noverlap,nfft,fs_egg);
                power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

                range_idx = find(w1>0.07 & w1<0.15); %Tachygastric range
                tachy_power(tasks,subs-1000,visit) = sum(power_norm(range_idx,:))./sum(power_norm); %Normogastric range

            end

            %Check whether file exists or not for bradygastric
            if ~isfile(append(brady_path,'V' , num2str(visit),'/',num2str(subs), ...
                    '/Task',num2str(tasks),'.mat'))
                fprintf('%d .mat doesnt exist\n',subs)
                brady_power(tasks,subs-1000,visit) = nan;
            else
                %Loading the egg signal
                load(append(brady_path,'V' , num2str(visit),'/',num2str(subs), ...
                    '/Task',num2str(tasks),'.mat'));
                %Zero padding the signal if the length is less than the
                %window size
                if length(egg_section) < 60*fs_egg
                    egg_section(length(egg_section):60*fs_egg) = 0; 
                end
                %PSD estimation 
                [Pxx1,w1] =  pwelch(egg_section,window,noverlap,nfft,fs_egg);
                power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

                range_idx = find(w1>0.0083 & w1<0.03); %Normogastric range
                brady_power(tasks,subs-1000,visit) = sum(power_norm(range_idx,:))./sum(power_norm); %Normogastric range

            end
        end
    end
end


%% Saving the variables
%save('EGG_power_v4.mat','normo_power','brady_power','tachy_power');       