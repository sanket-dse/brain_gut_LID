main_path_V1 = 'C:\Users\sanke\Downloads\EEG EXCEL DATA\V1\';
main_path_V3 = 'C:\Users\sanke\Downloads\EEG EXCEL DATA\V3\';
main_path = 'C:\Users\sanke\Downloads\PD_EEG_EGG\';
%subs_to_run = [1001:1028,1030:1065];
subs_to_run = 1001:1002;

for subs = subs_to_run 

    fileList_V1 = dir(append(main_path_V1 , num2str(subs)) + "\*.xlsx");
    fileList_V3 = dir(append(main_path_V3 , num2str(subs)) + "\*.xlsx");

    mkdir(append('C:\Users\sanke\Downloads\PD_EEG_EGG\mat_files\V1\',num2str(subs)));
    mkdir(append('C:\Users\sanke\Downloads\PD_EEG_EGG\mat_files\V3\',num2str(subs)));

    for files = 1:size(fileList_V1)
        
        path1 = append(main_path_V1,num2str(subs)) + "\" + fileList_V1(files,1).name;
        path2 = append(main_path_V3,num2str(subs)) + "\" + fileList_V3(files,1).name;
        
        [~,sheets1] = xlsfinfo(path1);
        [~,sheets2] = xlsfinfo(path2);
        
        data_V1=[];
        data_V3=[];
        for i=size(sheets1,2):-1:1
            fprintf('Loading sheet %d\n',size(sheets1,2)-(i-1));
            [numbers, strings, raw] = xlsread(path1,sheets1{:,i});
            data_V1 = cat(2,data_V1,transpose(numbers));
        end

        for i=size(sheets2,2):-1:1
            fprintf('Loading sheet %d\n',size(sheets2,2)-(i-1));
            [numbers, strings, raw] = xlsread(path2,sheets2{:,i});
            data_V3 = cat(2,data_V3,transpose(numbers));
        end

        save(main_path + "\mat_files\V1\" + fileList_V1(files,1).name(1:end-5) + ".mat",'data_V1');
        save(main_path + "\mat_files\V3\" + fileList_V3(files,1).name(1:end-5) + ".mat",'data_V3');

        
    end
end



