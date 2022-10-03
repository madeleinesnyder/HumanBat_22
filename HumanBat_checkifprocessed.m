function [processed_status] = HumanBat_checkifprocessed(batdate,batid)

% Function to check what parts of a dataset are processed
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',num2str(batid),'/processed/',num2str(batdate),'/');

% Cortex status
cortex_files = dir(fullfile(exp_data_path,'cortex/*.mat'));
for i=1:length(cortex_files)
    disp(cortex_files(i).name)
end
if exist(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'))
    disp("CORTEX FULLY PROCESSED per July 2022 pipeline stage");
end

% Ciholas status
ciholas_files = dir(fullfile(exp_data_path,'ciholas/*.mat'));
for i=1:length(ciholas_files)
    disp(ciholas_files(i).name)
end
if exist(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'))
    disp("CIHOLAS FULLY PROCESSED per July 2022 pipeline stage");
end

% Ephys status
for i=[13,15]
    logger_path = strcat(exp_data_path,'/ephys/logger',num2str(i),'/');
    single_units = dir(fullfile(logger_path,'/extracted_data/SingleUnits_*'));
    if ~isempty(single_units)
        disp(strcat("Logger ",num2str(i), " has SpikeSort3D sorted units."));
    elseif exist(strcat(logger_path,'/params.py'))
        disp(strcat("Logger ",num2str(i)," has been Kilosorted."));
    else
        disp(strcat("Logger ",num2str(i)," needs ephys sorting."));
    end
end

end