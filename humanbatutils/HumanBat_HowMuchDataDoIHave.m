function [data_struct] = HumanBat_HowMuchDataDoIHave(dates,loggers);

% Function for how much data I have on a given stretch of days and loggers

root_path = '/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/';

% Get number of cells on those days for those loggers
for i=1:length(dates)
    for j=1:length(loggers)
        cell_path = dir(fullfile(strcat(root_path,num2str(dates(i)),'/ephys/logger',num2str(loggers(j)),'/extracted_data/SingleUnit*')));
        load(strcat(cell_path.folder,'/',cell_path.name,'/','SingleUnits_',num2str(dates(i)),'.mat'));
        cell_mat(i,j) = length(TT_unit);
        clear TT_unit;
    end
end
HowMuchDataStruct.cells = cell_mat;
HowMuchDataStruct.cell_count = sum(cell_mat);

% Get number of flights on those days for each bat in the room
for i=1:length(dates)
    load(strcat(root_path,num2str(dates(i)),'/cortex/clustered_cortex_flights.mat'));
    flight_mat(i,1) = cortex_flights.N; 
    flight_clus(i,1) = size(cortex_flights.clusterIndex,2);
    clear cortex_flights;
    load(strcat(root_path,num2str(dates(i)),'/ciholas/clustered_ciholas_flights.mat'));
    flight_mat(i,2) = size(ciholas_flights.flights,1); 
    flight_clus(i,2) = length(unique(ciholas_flights.cluster_idx));
    clear ciholas_flights;
end
HowMuchDataStruct.cortexFlights = flight_mat(:,1);
HowMuchDataStruct.NcortexClusters = flight_clus(:,1);
HowMuchDataStruct.ciholasFlights = flight_mat(:,2);
HowMuchDataStruct.NciholasClusters = flight_clus(:,2);

% Get number of rewards
disp("Yikes. Todo");


end
