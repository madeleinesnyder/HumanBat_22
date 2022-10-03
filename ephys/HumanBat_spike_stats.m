function [] = HumanBat_spike_stats(dates,loggers)

% For a given set of date and loggers, compute how many putative units
% there are

server_root='/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/';

unit_count = 0;
for i=1:length(loggers)
    for j=1:length(dates)
        single_units_dirs=dir(fullfile(strcat(server_root,num2str(dates(j)),'/ephys/logger',num2str(loggers(i)),'/extracted_data/SingleUnits*')))
        load(strcat(single_units_dirs.folder,'/',single_units_dirs.name,'/SingleUnits_',num2str(dates(j)),'.mat'));
        unit_count = unit_count+length(TT_unit);
    end
end

disp(strcat("For logger(s) ",num2str(loggers)," on dates ",num2str(dates)," there are ",num2str(unit_count)," putative units"));

end