% Script for looking at cortex_flight_struct

% cortex_flight_struct should be loaded 
clear all;
batdate=220418;  logger=15;
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));

% Load the pre-pruned data if it exists
if exist(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_bat_final_flight_structure.mat'))
    disp("Loading pruned, resorted, and trimmed data");
    load(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_bat_final_flight_structure.mat'));
    data_type = 'prt';
elseif exist(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'))
    disp("Loading pruned and resorted data");
    load(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'));
    data_type = 'pr';
elseif exist(strcat(exp_data_path,'cortex/pruned_cortex_bat_final_flight_structure.mat'))
    disp("Loading pruned non-resorted data");
    load(strcat(exp_data_path,'cortex/pruned_cortex_bat_final_flight_structure.mat'));
    data_type = 'p';
else
    disp("Loading non-pruned, non-resorted data, CORTEX");
    load(strcat(exp_data_path,'cortex/cortex_bat_final_flight_structure.mat'));
    data_type='u';
end

%% === Do the flights look ok? Prune out the shitty ones for now
if strcmp(data_type,'u')
    for i=1:length(cortex_flight_struct)
        if isempty(cortex_flight_struct{i})
            continue;
        else
            figure(); hold on; scatter3(cortex_flight_struct{i}.pos(1,:),...
                cortex_flight_struct{i}.pos(2,:),...
                cortex_flight_struct{i}.pos(3,:));
            xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]); view(60,10); hold off;
            prompt = "Keep (K) or Toss (T)?";
            x = input(prompt,"s");
            if x=='T'
                cortex_flight_struct{i} = [];
                close all;
            else
                disp("Good flight, look at raster");
                close all;
            end
        end
    end
    save(strcat(exp_data_path,'cortex/pruned_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct','-v7.3');
else
    disp("Flights have already been pruned and sorted!")
end

%% === Create category matrix for all categories you want to raster
% Categories:
    % 2. All Flights To Madeleine = {indexes of flights in ciholas_flight_struct};
    % 3. All Flights To Kevin = ...
    % 1. All Flights Of Cluster X = {} {} {} ...

% Change variable names to get rid of the big cortex_flight_struct if it
% existed
to_resort=0;
if strcmp(data_type,'u') | strcmp(data_type,'p')
    to_resort = 1;
end

if strcmp(data_type,'u') | strcmp(data_type,'p')
    % 1. Clusters
    fclus_list = []; fclus_flight_idx = [];
    for i=1:length(cortex_flight_struct)
        if ~isempty(cortex_flight_struct{i})
            fclus_list = [fclus_list,cortex_flight_struct{i}.fclus];
            fclus_flight_idx = [fclus_flight_idx,i];
        end
    end
    clusters = unique(fclus_list); Flight_Group_Matrix_Clusters=cell(1,length(clusters));
    for j=1:length(clusters)
        for i=1:length(fclus_list)
            if fclus_list(i) == clusters(j)
                Flight_Group_Matrix_Clusters{j} = [Flight_Group_Matrix_Clusters{j},fclus_flight_idx(i)];
            end
        end
    end
    Flight_Group_Matrix.Clusters = Flight_Group_Matrix_Clusters; clear Flight_Group_Matrix_Clusters;

    % Plot a few clusters to check
    clus = 3;
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for i=1:length(Flight_Group_Matrix.Clusters{clus})
        flightnum = Flight_Group_Matrix.Clusters{clus}(i);
        plot3(cortex_flight_struct{flightnum}.pos(1,:),cortex_flight_struct{flightnum}.pos(2,:),cortex_flight_struct{flightnum}.pos(3,:));
    end
    title(strcat("Cluster ",num2str(clus)));
    hold off;
else
    % 1. Clusters
    fclus_list = []; fclus_flight_idx = [];
    for i=1:length(cortex_flight_struct_resort)
        if ~isempty(cortex_flight_struct_resort{i})
            fclus_list = [fclus_list,cortex_flight_struct_resort{i}.fclus];
            fclus_flight_idx = [fclus_flight_idx,i];
        end
    end
    clusters = unique(fclus_list); Flight_Group_Matrix_Clusters=cell(1,max(clusters));
    for j=1:length(clusters)
        cc = clusters(j);
        for i=1:length(fclus_list)
            if fclus_list(i) == clusters(j)
                Flight_Group_Matrix_Clusters{cc} = [Flight_Group_Matrix_Clusters{cc},fclus_flight_idx(i)];
            end
        end
    end
    Flight_Group_Matrix.Clusters = Flight_Group_Matrix_Clusters; clear Flight_Group_Matrix_Clusters;

    % Plot a few clusters to check
    clus = 3;
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for i=1:length(Flight_Group_Matrix.Clusters{clus})
        flightnum = Flight_Group_Matrix.Clusters{clus}(i);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:));
    end
    title(strcat("Cluster ",num2str(clus)));
    hold off;
end

% Resort the clusters
if to_resort == 1
    original_index = [1:length(cortex_flight_struct)];
    cluster_array = cell2mat(Flight_Group_Matrix.Clusters);
    og_idx = [];
    for i=1:length(original_index)
        og_idx_ = find(cluster_array==original_index(i));
        if isempty(og_idx_)
            og_idx{i} = NaN;
        else
            og_idx{i} = og_idx_;
        end
    end
end

% Get the list in order of cluster identities
cv = [];
for i=1:length(Flight_Group_Matrix.Clusters)
    cv_ = ones(1,length(Flight_Group_Matrix.Clusters{i}))*i;
    cv = [cv,cv_];
end

% Re-sort the clusters
for i=1:length(Flight_Group_Matrix.Clusters)
    Flight_Group_Matrix.Clusters{i} = find(cv==i);
end

if to_resort == 1
    cortex_flight_struct_resort = cortex_flight_struct(cluster_array');
    save(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');
end

% Get start times, end times CHRONOLOGICAL 
% Load in chronological struct (cortex_flight_struct)
%load(strcat(exp_data_path,'cortex/cortex_bat_final_flight_structure.mat'));

% if to_resort==1
%     chronological_start_times = []; chronological_end_times = [];
%     for i=1:length(cortex_flight_struct)
%         if ~isempty(cortex_flight_struct{i})
%             chronological_start_times = [chronological_start_times,cortex_flight_struct{i}.fstart/120];
%             chronological_end_times = [chronological_end_times,cortex_flight_struct{i}.fend/120];
%         end
%     end
% end

% 2. Madeleine/Kevin Flights     (cortex_flight_struct_resort)
cutoff = 700; cutoff_bat = 500; Flight_Group_Matrix_To_M=[]; Flight_Group_Matrix_To_K=[]; Flight_Group_Matrix_To_CiholasBat=[]; Flight_Group_Matrix_To_SelfSpot=[];
for i=1:length(cortex_flight_struct_resort)
    if ~isempty(cortex_flight_struct_resort{i})
        % Find if flight was to Madeleine
        mean_end_position = mean(cortex_flight_struct_resort{i}.pos(:,end-50:end)');
        mean_madeleine_position = mean(cortex_flight_struct_resort{i}.pos_M(end-50:end,:));
        mean_kevin_position = mean(cortex_flight_struct_resort{i}.pos_K(end-50:end,:));
        try
            mean_CiholasBat_position = mean(cortex_flight_struct_resort{i}.pos_ciholas_bat(end-20:end,:));
            ciholas_bat_neg=0;
        catch
            disp(strcat("No CiholasBat Position on flight ",num2str(i)));
            ciholas_bat_neg=1;
        end
%             figure(); hold on; 
%             scatter3(cortex_flight_struct{i}.pos(:,1),cortex_flight_struct{i}.pos(:,2),cortex_flight_struct{i}.pos(:,3));
%             scatter3(mean_madeleine_position(1),mean_madeleine_position(2),mean_madeleine_position(3),60); 
%             scatter3(mean_kevin_position(1),mean_kevin_position(2),mean_kevin_position(3),60); 
%             scatter3(mean_CortexBat_position(1),mean_CortexBat_position(2),mean_CiholasBat_position(3)); 
%             scatter3(mean_end_position(1),mean_end_position(2),mean_end_position(3),100,'filled');
%             title('Flight ',num2str(i)); hold off;

        if pdist2(mean_end_position, mean_madeleine_position, 'euclidean') < cutoff
            Flight_Group_Matrix_To_M = [Flight_Group_Matrix_To_M,i];
        elseif pdist2(mean_end_position, mean_kevin_position, 'euclidean') < cutoff
            Flight_Group_Matrix_To_K = [Flight_Group_Matrix_To_K,i];
        elseif ciholas_bat_neg==0
            if pdist2(mean_end_position, mean_CiholasBat_position, 'euclidean') < cutoff_bat
                Flight_Group_Matrix_To_CiholasBat = [Flight_Group_Matrix_To_CiholasBat,i];
            end
        elseif ciholas_bat_neg==1
            Flight_Group_Matrix_To_CiholasBat = [Flight_Group_Matrix_To_CiholasBat];
        else
            Flight_Group_Matrix_To_SelfSpot = [Flight_Group_Matrix_To_SelfSpot,i];
        end
    end
end

Flight_Group_Matrix.AllFlightsToMadeleine = Flight_Group_Matrix_To_M; clear Flight_Group_Matrix_To_M;
Flight_Group_Matrix.AllFlightsToKevin = Flight_Group_Matrix_To_K; clear Flight_Group_Matrix_To_K;
Flight_Group_Matrix.AllFlightsToCiholasBat = Flight_Group_Matrix_To_CiholasBat; clear Flight_Group_Matrix_To_CiholasBat;
Flight_Group_Matrix.Flight_Group_Matrix_To_SelfSpot = Flight_Group_Matrix_To_SelfSpot; clear Flight_Group_Matrix_To_SelfSpot;

% load in extra data needed
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
%load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));

% Get start times, end times CLUSTER SORTED 
cluster_sorted_start_times = []; cluster_sorted_end_times = [];
for i=1:length(cortex_flight_struct_resort)
    if ~isempty(cortex_flight_struct_resort{i})
        cluster_sorted_start_times = [cluster_sorted_start_times,cortex_flight_struct_resort{i}.fstart/120];
        cluster_sorted_end_times = [cluster_sorted_end_times,cortex_flight_struct_resort{i}.fend/120];
        %cluster_sorted_start_times = [cluster_sorted_start_times,cortex_flight_struct_resort{i}.fstart_trimmed/120];
        %cluster_sorted_end_times = [cluster_sorted_end_times,cortex_flight_struct_resort{i}.fend_trimmed/120];
    end
end

% Get cv
cv = [];
for i=1:length(cortex_flight_struct_resort)
    cv = [cv,cortex_flight_struct_resort{i}.fclus];
end

% Specify if trimmed or not
if strcmp(data_type,'prt')
    trimmed=1;
else
    trimmed=0;
end

    %% Raster plots of all kinds

    % -----------------------------------------------------------
    %                  Aligned to Takeoff            Aligned to Landing
    % Chronological
    % By Cluster
    % From Human 
    % To Human
    % From Location
    % To Location
    
    %% Chronological; Takeoff
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        % Take all usec times of flight takeoff
        t_temp_takeoff = chronological_start_times;
        
        % Add spices
        g_id = ones(1,length(t_temp_takeoff));
        clr = repmat(lines(1),length(t_temp_takeoff),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s',t_temp_takeoff',g_id,clr,interval_,[],strcat("Spikes from Cortex Bat Unit ",num2str(ii)," aligned to takeoff of Cortex flight CHRONOLOGICAL"),1,1);
    end
    
    %% Chronological; Landing
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        % Take all usec times of flight takeoff
        t_temp_takeoff = chronological_end_times;
        
        % Add spices
        g_id = ones(1,length(t_temp_takeoff));
        clr = repmat(lines(1),length(t_temp_takeoff),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s',t_temp_takeoff',g_id,clr,interval_,[],strcat("Spikes from Cortex Bat Unit ",num2str(ii)," aligned to landing of Cortex flight CHRONOLOGICAL"),1,1);
    end
    
    %% By Cluster; Takeoff
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        % Take all usec times of flight takeoff
        t_temp_takeoff = cluster_sorted_start_times;
        
        % Add spices
        g_id = ones(1,length(t_temp_takeoff));
        clr = repmat(lines(1),length(t_temp_takeoff),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s',t_temp_takeoff',g_id,clr,interval_,cv,strcat("Spikes from Cortex Bat Unit ",num2str(ii)," aligned to takeoff of Cortex flights sorted by cluster identity"),1,1);
    end
    
    %% By Cluster; Landing 
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        % Take all usec times of flight takeoff
        t_temp_takeoff = cluster_sorted_end_times;
        
        % Add spices
        g_id = ones(1,length(t_temp_takeoff));
        clr = repmat(lines(1),length(t_temp_takeoff),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s',t_temp_takeoff',g_id,clr,interval_,cv,strcat("Spikes from Cortex Bat Unit ",num2str(ii)," aligned to landing of Cortex flights sorted by cluster identity"),1,1);
    end
    
    %% Plot all flights to Madeleine and all flights to Kevin 
    figure(); 
    subplot(1,2,1); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for k=1:length(Flight_Group_Matrix.AllFlightsToKevin)
        flightnum = Flight_Group_Matrix.AllFlightsToKevin(k);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:));
    end
    title("All Flights to Kevin");
    hold off;
    subplot(1,2,2); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for k=1:length(Flight_Group_Matrix.AllFlightsToMadeleine)
        flightnum = Flight_Group_Matrix.AllFlightsToMadeleine(k);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:));
    end
    title("All Flights to Madeleine");
    hold off; 
    
    %% To Madeleine versus To Kevin (sorted by cluster); Takeoff 
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        to_list_m = []; to_list_k = [];
        for j=1:length(Flight_Group_Matrix.AllFlightsToMadeleine)
            idx = Flight_Group_Matrix.AllFlightsToMadeleine(j);
            to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart/120];
            %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
        end
        t_temp_takeoff_m = to_list_m'; clear to_list_m;
        cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
        for j=1:length(Flight_Group_Matrix.AllFlightsToKevin)
            idx = Flight_Group_Matrix.AllFlightsToKevin(j);
            to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart/120];
            %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
        end
        t_temp_takeoff_k = to_list_k'; clear to_list_k;
        cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);
    
        % Add spices
        g_id_m = ones(1,length(t_temp_takeoff_m));     g_id_k = ones(1,length(t_temp_takeoff_k));
        clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);     clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,cv_m,cv_k,strcat("CiholasBatUnit ",num2str(ii)," : takeoff of cortexFlights To M"),strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To K"),1,1);
    end
    
    %% To Madeleine versus to Kevin (sorted by cluster); Landing
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
       to_list_m = []; to_list_k = [];
        for j=1:length(Flight_Group_Matrix.AllFlightsToMadeleine)
            idx = Flight_Group_Matrix.AllFlightsToMadeleine(j);
            to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend/120];
            %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
        end
        t_temp_takeoff_m = to_list_m'; clear to_list_m;
        cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
        for j=1:length(Flight_Group_Matrix.AllFlightsToKevin)
            idx = Flight_Group_Matrix.AllFlightsToKevin(j);
            %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend/120];
        end
        t_temp_takeoff_k = to_list_k'; clear to_list_k;
        cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);
    
        % Add spices
        g_id_m = ones(1,length(t_temp_takeoff_m));     g_id_k = ones(1,length(t_temp_takeoff_k));
        clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);     clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,cv_m,cv_k,strcat("cortexBatUnit ",num2str(ii)," : landing of cortexFlights To M"),strcat("cortexBatUnit ",num2str(ii)," : landing of cortexFlights To K"),1,1);
    end
    
    %% === Cluster X Human === %%
    Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
    Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
    Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
    Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)
    
    %% To Madeleine vs Kevin x Cluster (clus); Takeoff 
    for i=1:length(sig_units_dur_clus{clus})%length(B_ephys_data.TT_unit)
        ii= sig_units_dur_clus{clus}(i); %unit;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        % Take all usec times of flight takeoff CHANGE THIS FOR THE
        % CATEGORIES!!!!
        to_list_m = [];
        for j=1:length(Mx_flights)
            idx = Mx_flights(j);
            %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart/120];
        end
        t_temp_takeoff_m = to_list_m'; clear to_list;
        to_list_k = [];
        for j=1:length(Kx_flights)
            idx = Kx_flights(j);
            %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart/120];
        end
        t_temp_takeoff_k = to_list_k'; clear to_list;
        cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
        cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);
    
        % Add spices
        g_id_m = ones(1,length(t_temp_takeoff_m));
        clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);
        g_id_k = ones(1,length(t_temp_takeoff_k));
        clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To M"),strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To K"),1,1);
    end
    
    %% To Madeleine v.s Kevin x Cluster (clus); Landing
    for i=1:length(sig_units_post_clus{clus})%length(B_ephys_data.TT_unit)
        ii=sig_units_post_clus{clus}(i)%i%unit;%cool_list(i);
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
         % Take all usec times of flight takeoff CHANGE THIS FOR THE
        % CATEGORIES!!!!
        to_list_m = [];
        for j=1:length(Mx_flights)
            idx = Mx_flights(j);
            %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend/120];
        end
        t_temp_takeoff_m = to_list_m'; clear to_list;
        to_list_k = [];
        for j=1:length(Kx_flights)
            idx = Kx_flights(j);
            %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend/120];
        end
        t_temp_takeoff_k = to_list_k'; clear to_list;
        cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
        cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);
    
        % Add spices
        g_id_m = ones(1,length(t_temp_takeoff_m));
        clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);
        g_id_k = ones(1,length(t_temp_takeoff_k));
        clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("cortexBatUnit ",num2str(ii)," : landing of cortexFlights To M"),strcat("cortexBatUnit ",num2str(ii)," : landing of cortexFlights To K"),1,1);
    end
    
    %% Plot in-cluster flights to Kevin versus Madeleine 
    figure(); hold on;
    for k=1:length(Flight_Group_Matrix.Clusters{clus})
        flightnum = Flight_Group_Matrix.Clusters{clus}(k);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:),'r');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,1),cortex_flight_struct_resort{flightnum}.pos(2,1),cortex_flight_struct_resort{flightnum}.pos(3,1),50,'g','filled');     
    end
    
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for k=1:length(Kx_flights)
        flightnum = Kx_flights(k);%og_idx{Kx_flights(k)};
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:),'r');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,1),cortex_flight_struct_resort{flightnum}.pos(2,1),cortex_flight_struct_resort{flightnum}.pos(3,1),50,'g','filled'); 
    end
    for k=1:length(Mx_flights)
        flightnum = Mx_flights(k);%og_idx{Mx_flights(k)};
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:),'b');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,1),cortex_flight_struct_resort{flightnum}.pos(2,1),cortex_flight_struct_resort{flightnum}.pos(3,1),50,'g','filled'); 
    end
    title(strcat("All Cluster ",num2str(clus)," Flights to K (red) Flights to M (blue)"));
    hold off; 

    % Plot which rasters are on which flights 
    HumanBat_joint_raster_flightpath_plots(cortex_flight_struct_resort,B_ephys_data,Mx_flights,Kx_flights)

    %% To Location x versus from Location x 
    flights_to_smush = '_8'; % (From HumanBat_SI.mat)
    
    % If it is a "to" flight...
    Location_Flight_List = [];
    if flights_to_smush(1) == '_'
        for j=1:10
            if str2double(flights_to_smush(2)) == j
                index_ = j
                disp(strcat("Looking at rasters of takeoffs of all flights TO location ",num2str(index_)," Madeleine versus Kevin"));
            end
        end 
        % Find all the flights in this Location subset
        for i=1:size(cortex_flight_struct_resort,2)
            if cortex_flight_struct_resort{i}.tripod_landing == index_
                Location_Flight_List = [Location_Flight_List,i];
            end
        end
    elseif flights_to_smush(2) == '_'
        for j=1:10
            if str2double(flights_to_smush(1)) == j
                index_ = j
                disp(strcat("Looking at rasters of takeoffs of all flights FROM location ",num2str(index_)," Madeleine versus Kevin"));
            end
        end
        % Find all the flights in this Location subset
        for i=1:size(cortex_flight_struct_resort,2)
            if cortex_flight_struct_resort{i}.tripod_takeoff == index_
                Location_Flight_List = [Location_Flight_List,i];
            end
        end
    end
    Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Location_Flight_List);
    Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
    Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Location_Flight_List);
    Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)
    Void_flights = Location_Flight_List;

    %% To/From Location X. Madeleine versus To Kevin; Takeoff 
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        to_list_m = []; to_list_k = [];
        for j=1:length(Mx_flights)
            idx = Mx_flights(j);
            %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart/120];
        end
        t_temp_takeoff_m = to_list_m'; clear to_list_m;
        cv_m = cv(Mx_flights);
        for j=1:length(Kx_flights)
            idx = Kx_flights(j);
            %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart/120];
        end
        t_temp_takeoff_k = to_list_k'; clear to_list_k;
        cv_k = cv(Kx_flights);
    
        % Add spices
        g_id_m = ones(1,length(t_temp_takeoff_m));     g_id_k = ones(1,length(t_temp_takeoff_k));
        clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);     clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,cv_m,cv_k,strcat("CiholasBatUnit ",num2str(ii)," : takeoff of ",flights_to_smush," CiholasFlights To M"),strcat("CiholasBatUnit ",num2str(ii)," : takeoff of ",flights_to_smush," CiholasFlights To K"),1,1);
    end
    
    %% To/From Location X. Madeleine versus To Kevin; Landing 
    for i=1:length(B_ephys_data.TT_unit)
        ii=i;
        s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
        
        to_list_m = []; to_list_k = [];
        for j=1:length(Mx_flights)
            idx = Mx_flights(j);
            %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend/120];
        end
        t_temp_takeoff_m = to_list_m'; clear to_list_m;
        cv_m = cv(Mx_flights);
        for j=1:length(Kx_flights)
            idx = Kx_flights(j);
            %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend/120];
        end
        t_temp_takeoff_k = to_list_k'; clear to_list_k;
        cv_k = cv(Kx_flights);
    
        % Add spices
        g_id_m = ones(1,length(t_temp_takeoff_m));     g_id_k = ones(1,length(t_temp_takeoff_k));
        clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);     clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
        interval_ = [-3 5];
    
        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,cv_m,cv_k,strcat("CiholasBatUnit ",num2str(ii)," : takeoff of ",flights_to_smush," CortexFlights To M"),strcat("CortexBatUnit ",num2str(ii)," : landing of ",flights_to_smush," CortexFlights To K"),1,1);
    end
    
    %% Plot x Location flights to Kevin versus Madeleine 
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for k=1:length(Void_flights)
        flightnum = Void_flights(k);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:),'Color',[0.6 0.6 0.6],'LineStyle',':');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,1),cortex_flight_struct_resort{flightnum}.pos(2,1),cortex_flight_struct_resort{flightnum}.pos(3,1),50,'g','filled');    
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,end),cortex_flight_struct_resort{flightnum}.pos(2,end),cortex_flight_struct_resort{flightnum}.pos(3,end),50,'m','filled');
    end
    for k=1:length(Kx_flights)
        flightnum = Kx_flights(k);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:),'r');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,1),cortex_flight_struct_resort{flightnum}.pos(2,1),cortex_flight_struct_resort{flightnum}.pos(3,1),50,'g','filled');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,end),cortex_flight_struct_resort{flightnum}.pos(2,end),cortex_flight_struct_resort{flightnum}.pos(3,end),50,'m','filled');
    end
    for k=1:length(Mx_flights)
        flightnum = Mx_flights(k);
        plot3(cortex_flight_struct_resort{flightnum}.pos(1,:),cortex_flight_struct_resort{flightnum}.pos(2,:),cortex_flight_struct_resort{flightnum}.pos(3,:),'b');
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,1),cortex_flight_struct_resort{flightnum}.pos(2,1),cortex_flight_struct_resort{flightnum}.pos(3,1),50,'g','filled');    
        scatter3(cortex_flight_struct_resort{flightnum}.pos(1,end),cortex_flight_struct_resort{flightnum}.pos(2,end),cortex_flight_struct_resort{flightnum}.pos(3,end),50,'m','filled');
    end
    title(strcat("All ",flights_to_smush," flights: K (red) Flights to M (blue)"));
    hold off;  

