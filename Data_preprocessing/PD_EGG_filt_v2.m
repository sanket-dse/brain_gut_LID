addpath(genpath('/data3/Sanket/eeglab2023.1'));
egg_path = '/data4/NEURO_DrRathour/DATA/OPEN BCI/';
filt_path = '/data3/Sanket/';
visits = [1];
%subs_to_run = [1001,1003:1028,1030:1049,1051:1067]; %subs 1002,1029 and 1050 not being run 
subs_to_run = [1022]; %subs 1002,1029 and 1050 not being run 
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
       
        %{ 
        % Filter signal
        fprintf('Being filtered\n')
        s1_filt = bandpass(egg_sig,[0.03 0.07],250);
        fprintf('FIltering done\n')
        %}

        % prepare the filter: finite impulse response filter with a bandwith of 
        % +/- 0.08 Hz of the peak frequency
        srate               = 250;
        center_frequency    = 0.07;        
        bandwidth           = 0.08;
        transition_width    = 0.15;
        nyquist             = srate/2;
        ffreq(1)            = 0;
        ffreq(2)            = (1-transition_width)*(bandwidth - center_frequency);
        ffreq(3)            = (bandwidth - center_frequency);
        ffreq(4)            = (center_frequency+bandwidth);
        ffreq(5)            = (1+transition_width)*(center_frequency+bandwidth);
        ffreq(6)            = nyquist;
        ffreq               = ffreq/nyquist;
        fOrder              = 3; % in cycles
        filterOrder         = fOrder*fix(srate/(bandwidth - center_frequency)); %in samples
        idealresponse       = [ 0 0 1 1 0 0 ];
        filterweights       = fir2(filterOrder,ffreq,idealresponse);

        % filter
        %EGG_filt = EGG_raw;
        egg_length = size(egg_sig,1);
        egg_sig_padded = [egg_sig ; zeros(225100 - egg_length,4)];
        disp('Filtering EGG - this will take some time');
        s1_filt   = filtfilt(filterweights,1,egg_sig_padded);
        fprintf('filtered');
        s1_filt = s1_filt(1:egg_length,:);

        %Create a directory
        mkdir(append(filt_path,'filtered_EGG_v2/V',num2str(visit),'/'));
        save(append(filt_path,'filtered_EGG_v2/V',num2str(visit),'/',num2str(sub),'.mat'),'s1_filt'); %saving the file
    end
end