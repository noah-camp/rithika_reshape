clear all
close all

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
warning('off', 'MATLAB:table:RowsAddedExistingVars');

data_path = 'D:\SleepDB\Intermediate';
oldcsv_path = 'D:\SleepDB\Intermediate_csv';
newcsv_path = 'D:\SleepDB\Output_csvs';
output_path = 'D:\SleepDB\Output_edfs';
mat_path = 'D:\SleepDB\Mats';
edfs = dir(fullfile(data_path, '*.edf'));
csvs = dir(fullfile(oldcsv_path, '*.csv'));
mkdir(output_path);
mkdir(newcsv_path);
mkdir(mat_path);

for i = 11:length(edfs)
    fprintf('Converting %i out of %i files (%s)...', i, length(edfs), edfs(i).name);
    [header signalHeader signalCell annotations] = blockEdfLoadAnno(strcat(data_path,'\', edfs(i).name));
    
    if(length(signalCell) ~= 53)
        fileID = fopen('irregular_files.txt','w');
        fprintf(fileID, '%s\n', edfs(i).name); 
        fclose(fileID);
        warning('File does not have 53 channels.');
        continue;
    end
    
    mkdir(strcat(edfs(i).folder, '\',(edfs(i).name(1:end-4)), '_EDF'))
    
    signalFlags = [false false];
    for j = 1:length(signalHeader)
        signalHeader(j).samples_in_record = signalHeader(j).samples_in_record/2;
        if(strcmp(signalHeader(j).signal_labels, 'LAT1') || strcmp(signalHeader(j).signal_labels, 'ECG'))
            signalHeader(j).signal_labels = 'ECG';
            ecg_chan = j;
            signalFlags(1) = true;
        elseif(strcmp(signalHeader(j).signal_labels, 'SpO2'))
            signalHeader(j).signal_labels ='Not SpO2';
        elseif(strcmp(signalHeader(j).signal_labels, 'OSAT') || strcmp(signalHeader(j).signal_labels, 'SpO2 B-B'))
            signalHeader(j).signal_labels ='SpO2';
            spo2_chan = j;
            signalFlags(2) = true;
        end
    end
    
    if(~signalFlags(1))
        warning('ECG Signal not found for %s', edfs(i).name );
    end
    if(~signalFlags(2))
        warning('SpO2 Signal not found for %s', edfs(i).name );
    end
    
    if(header.data_record_duration ~= 1)
        header.data_record_duration = 1;
        header.num_data_records = header.num_data_records * 2;
    end
    
    % Reduce sampling rate of SpO2
    orig_spo2 = signalCell{51};
    orig_spo2_sr = signalHeader(spo2_chan).samples_in_record;
    
    if(signalHeader(spo2_chan).samples_in_record ~= 1)
        signalCell{spo2_chan} = signalCell{spo2_chan}(1:signalHeader(spo2_chan).samples_in_record:end);
        signalHeader(spo2_chan).samples_in_record = 1;
    end

    status = blockEdfWrite(strcat(edfs(i).folder, '\',(edfs(i).name(1:end-4))), header, signalHeader, signalCell);
    %    movefile(strcat(data_path,'\', edfs(i).name), output_path);
    
    t = readtable(strcat(oldcsv_path,'\', csvs(i).name));
    sleepStages = t.DefaultStagingSet_stage_;
    %
    %    numOx = sum(strcmp(annotations(:, 2), "Oxygen Desaturation"));
    
    pptid = repmat("", size(sleepStages, 1)-1, 1);
    epoch = zeros(size(sleepStages, 1)-1, 1);
    EpochStartTime = zeros(size(sleepStages, 1)-1, 1);
    SleepStage = zeros(size(sleepStages, 1)-1, 1);
    Name = repmat("", size(sleepStages, 1)-1, 1);
    Start = repmat("", size(sleepStages, 1)-1, 1);
    Duration = repmat("", size(sleepStages, 1)-1, 1);
    Input = repmat("", size(sleepStages, 1)-1, 1);
    LowestSpO2 = repmat("", size(sleepStages, 1)-1, 1);
    Desaturation =repmat("", size(sleepStages, 1)-1, 1);
    T = table(pptid, epoch, EpochStartTime, SleepStage, Name, Start, Duration, Input, LowestSpO2, Desaturation);
%     T.pptid = repmat({''}, size(sleepStages, 1)-1, 1);
    
    
    % convert sleep stages
    for j=1:length(sleepStages)
        switch sleepStages(j)
            case 1
                sleepStages(j) = 3;
            case 2
                sleepStages(j) = 2;
            case 3
                sleepStages(j) = 1;
            case 4
                sleepStages(j) = 5;
            case 5
                sleepStages(j) = 0;
        end
    end
    
    % Create noise files for Dan's PPG software
    sleepStages_temp = sleepStages;
    sleepStages = sleepStages(1:end-1);
    ecg = signalCell{ecg_chan};
    ecgSamplingRate = signalHeader(ecg_chan).samples_in_record;
    ppg = orig_spo2;
    ppgSamplingRate = orig_spo2_sr;
    scoring_epoch_size_sec = 30;
    save(strcat(mat_path, '\', strrep(edfs(i).name, '.EDF', '.mat')), 'ecg', 'ecgSamplingRate', 'ppg', 'ppgSamplingRate', 'scoring_epoch_size_sec', 'sleepStages')
    sleepStages = sleepStages_temp;
    
    % setup variables
    keys = ["Oxygen Desaturation", "EEG arousal", "Obstructive Hypopnea"];
    replace = ["SpO2 artifact", "Arousal (ASDA)", "Hypopnea"];
    epoch = 0;
    EpochStartTime = -30;
    tablePos = 1;
    sleepPos = 2;
    flag = true;
    
    
    for j = 1:size(annotations, 1)
        % break out if no more sleep stages
        if(sleepPos > size(sleepStages, 1))
            break;
        end
        
        % break out if there's no more recording
        if EpochStartTime+30 >= length(signalCell{spo2_chan})
            break;
        end
        
        
        %----------------------------------- if key is found in annotations
        if contains(annotations(j, 2), keys)
            loc = annotations(j, 2) == keys;
            name = replace(loc);
            time_dur = double(strsplit(annotations(j, 1), char(21)));
            
            
            % ------------------- if this is a new epoch and a key is found
            if  abs(EpochStartTime - time_dur(1)) >= 30
                T.EpochStartTime(tablePos) = EpochStartTime + 30;
                EpochStartTime = T.EpochStartTime(tablePos);
                T.epoch(tablePos) = epoch + 1;
                epoch = T.epoch(tablePos);
                T.SleepStage(tablePos) = sleepStages(sleepPos);
                T.Name(tablePos) = name;
                T.Start(tablePos) = {num2str(round(time_dur(1), 1))};
                T.Duration(tablePos) = {num2str(round(time_dur(2), 1))};
                if loc(1)
                    T.Input(tablePos) = "SpO2";
                    start = round(time_dur(1));
                    stop = round((time_dur(1) + time_dur(2)));
                    minO2 = min(signalCell{spo2_chan}(start:stop));
                    maxO2 = max(signalCell{spo2_chan}(start:stop));
                    T.LowestSpO2(tablePos) = {num2str(round(minO2))};
                    T.Desaturation(tablePos) = {num2str(round(maxO2) - round(minO2))};
                elseif loc(2)
                    T.Input(tablePos) = "EEG1";
                    T.LowestSpO2(tablePos) = "NA";
                    T.Desaturation(tablePos) = "NA";
                else
                    T.Input(tablePos) = "NA";
                    T.LowestSpO2(tablePos) = "NA";
                    T.Desaturation(tablePos) = "NA";
                end
                
                sleepPos = sleepPos+1;
                flag = false;
            
                
            % ------------------- if this the same epoch and a key is found 
            else
                if(flag)
                    tablePos = tablePos-1;
                    flag = false;
                end
                
                T.EpochStartTime(tablePos) = EpochStartTime;
                T.epoch(tablePos) = epoch;
                T.SleepStage(tablePos) = sleepStages(sleepPos);
                T.Name(tablePos) = name;
                T.Start(tablePos) = {num2str(round(time_dur(1), 1))};
                T.Duration(tablePos) = {num2str(round(time_dur(2), 1))};
                if loc(1)
                    T.Input(tablePos) = "SpO2";
                    start = round(time_dur(1));
                    stop = round(time_dur(1) + time_dur(2));
                    minO2 = min(signalCell{spo2_chan}(start:stop));
                    maxO2 = max(signalCell{spo2_chan}(start:stop));
                    T.LowestSpO2(tablePos) = {num2str(round(minO2))};
                    T.Desaturation(tablePos) = {num2str(round(maxO2) - round(minO2))};
                elseif loc(2)
                    T.Input(tablePos) = "EEG1";
                    T.LowestSpO2(tablePos) = "NA";
                    T.Desaturation(tablePos) = "NA";
                else
                    T.Input(tablePos) = "NA";
                    T.LowestSpO2(tablePos) = "NA";
                    T.Desaturation(tablePos) = "NA";
                end
            end
            tablePos = tablePos+1;
        
            
        % ---------------------- if this is a new epoch and no key is found
        elseif abs(double(annotations(j, 1)) - EpochStartTime) >= 30
            
            T.EpochStartTime(tablePos) = EpochStartTime + 30;
            EpochStartTime = T.EpochStartTime(tablePos);
            T.epoch(tablePos) = epoch + 1;
            epoch = T.epoch(tablePos);
            T.SleepStage(tablePos) = sleepStages(sleepPos);
            T.Name(tablePos) = "NA";
            T.Start(tablePos) = "NA";
            T.Duration(tablePos) = "NA";
            T.Input(tablePos) = "NA";
            T.LowestSpO2(tablePos) = "NA";
            T.Desaturation(tablePos) = "NA";

            sleepPos = sleepPos+1;
            tablePos = tablePos+1;
            flag = true;
            
            w = warning('query','last');

        end
    end
    
    T(:, 1) = {edfs(i).name};    
    writetable(T,strrep(strcat(newcsv_path, '\', csvs(i).name), '_Hypnogram', '.edf'),'QuoteStrings',true)
    
    clear signalCell
    clear signalHeader
    clear header
end