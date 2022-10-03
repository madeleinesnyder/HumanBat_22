function [Unit_Summary] = HumanBat_UnitSummary_LandingLocations_Cortex(batdate,logger,ciholas_or_cortex,unit,cortex_flight_struct_resort,B_ephys_data)

% Function that makes a unit summary figure for a given unit
% Heatmap of Spatial Information
% Red dot map of Spatial information
% Spatial Information Calculation
% Which clusters had enough flights, plotted in red and blue
% Rasters for MvsK on those clusters in red and blue.
% Significance shuffle to show how out-of-distribution the difference
% between M and K are for that unit and cluster.

% Exp data path:
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');

% Load in the ephys data:
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));

% Load in the pruned, resorted, trimmed ciholas flight data:
load(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_bat_final_flight_structure.mat'));


% Load in the data path 
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');

% DEFINE: RDC_SPLIT (Whic cluster are you looking at?)
rdc_split_choice = 2;

% Make spatial information figures
if ciholas_or_cortex == 'cortex'

    % SI info and heatmaps
    SI_info = HumanBat_SI_fx_Cortex(batdate,logger,unit,'all');
    Unit_Summary.SI_info = SI_info;

    % Figure out how many flights per cluster
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

    red_dot_clusterable = [];
    for i=1:length(Flight_Group_Matrix.Clusters)
        if length(Flight_Group_Matrix.Clusters{i}) > 5 
            red_dot_clusterable = [red_dot_clusterable,i];
        end
    end

    % Red dot figures
    HumanBat_SI_red_dots_Cortex(batdate,logger,unit,cortex_flight_struct_resort,red_dot_clusterable);

    % Madeleine versus Kevin figures
    cutoff = 800; cutoff_bat = 500; Flight_Group_Matrix_To_M=[]; Flight_Group_Matrix_To_K=[]; Flight_Group_Matrix_To_CortexBat=[]; Flight_Group_Matrix_To_SelfSpot=[];
    for i=1:length(cortex_flight_struct_resort)
        if ~isempty(cortex_flight_struct_resort{i})
            % Find if flight was to Madeleine
            mean_end_position = mean(cortex_flight_struct_resort{i}.pos(:,end-50:end)');
            mean_madeleine_position = mean(cortex_flight_struct_resort{i}.pos_M(end-50:end,:));
            mean_kevin_position = mean(cortex_flight_struct_resort{i}.pos_K(end-50:end,:));
            mean_CortexBat_position = mean(cortex_flight_struct_resort{i}.pos_ciholas_bat(end-20:end,:));
    %             figure(); hold on; 
    %             scatter3(cortex_flight_struct{i}.pos(:,1),cortex_flight_struct{i}.pos(:,2),cortex_flight_struct{i}.pos(:,3));
    %             scatter3(mean_madeleine_position(1),mean_madeleine_position(2),mean_madeleine_position(3),60); 
    %             scatter3(mean_kevin_position(1),mean_kevin_position(2),mean_kevin_position(3),60); 
    %             scatter3(mean_CortexBat_position(1),mean_CortexBat_position(2),mean_CortexBat_position(3)); 
    %             scatter3(mean_end_position(1),mean_end_position(2),mean_end_position(3),100,'filled');
    %             title('Flight ',num2str(i)); hold off;
    
            if pdist2(mean_end_position, mean_madeleine_position, 'euclidean') < cutoff
                Flight_Group_Matrix_To_M = [Flight_Group_Matrix_To_M,i];
            elseif pdist2(mean_end_position, mean_kevin_position, 'euclidean') < cutoff
                Flight_Group_Matrix_To_K = [Flight_Group_Matrix_To_K,i];
            elseif pdist2(mean_end_position, mean_CortexBat_position, 'euclidean') < cutoff_bat
                Flight_Group_Matrix_To_CortexBat = [Flight_Group_Matrix_To_CortexBat,i];
            else
                Flight_Group_Matrix_To_SelfSpot = [Flight_Group_Matrix_To_SelfSpot,i];
            end
        end
    end
    
    Flight_Group_Matrix.AllFlightsToMadeleine = Flight_Group_Matrix_To_M; clear Flight_Group_Matrix_To_M;
    Flight_Group_Matrix.AllFlightsToKevin = Flight_Group_Matrix_To_K; clear Flight_Group_Matrix_To_K;
    Flight_Group_Matrix.AllFlightsToCortexBat = Flight_Group_Matrix_To_CortexBat; clear Flight_Group_Matrix_To_CortexBat;
    Flight_Group_Matrix.Flight_Group_Matrix_To_SelfSpot = Flight_Group_Matrix_To_SelfSpot; clear Flight_Group_Matrix_To_SelfSpot;

    % For each cluster that has enough flights, look at the M-K split
    rdc_split = [];
    for i=1:length(red_dot_clusterable)
        clus = red_dot_clusterable(i);
        Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
        Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
        Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
        Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)

        if length(Mx_flights) > 2 & length(Kx_flights) > 2
            disp(strcat("Cluster ",num2str(clus), " has enough flights to look at MvK split"));
            rdc_split = [rdc_split,clus];
            figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
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
            title(strcat("All Cluster ",num2str(clus)," Flights to K (red) Flights to M (blue)"));
            hold off; 
        end
    end
    Unit_Summary.rdc = rdc_split;

    % Get cv
    cv = [];
    for i=1:length(cortex_flight_struct_resort)
        cv = [cv,cortex_flight_struct_resort{i}.fclus];
    end

    % Figure out if this unit has sig diff activity for M versus K flights
    % If it does have sig pre, or sig landing, plot the rasters and cluster
    disp("HAVE NOT MODIFIED FOR CORTEX HEREON OUT");
    clear sig_pre sig_land
    for i=rdc_split_choice:rdc_split_choice%length(rdc_split)
        clus = rdc_split(i);
        Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
        Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
        Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
        Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)

        for j=1:5
        [sig_pre(j),sig_land(j)] = HumanBat_MvK_Permutation_Clusters_Cortex(clus,unit,Flight_Group_Matrix,cortex_flight_struct_resort,B_ephys_data)
        end

        if mean(sig_pre) > 0.89
            s = B_ephys_data.TT_unit(unit).AlignedTimestamps/1e6;
            
            % Take all usec times of flight takeoff CHANGE THIS FOR THE
            % CATEGORIES!!!!
            to_list_m = [];
            for j=1:length(Mx_flights)
                idx = Mx_flights(j);
                to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            end
            t_temp_takeoff_m = to_list_m'; clear to_list;
            to_list_k = [];
            for j=1:length(Kx_flights)
                idx = Kx_flights(j);
                to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
            end
            t_temp_takeoff_k = to_list_k'; clear to_list;
            cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
            cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);
        
            to_list_m = [];
            for j=1:length(Mx_flights)
                idx = Mx_flights(j);
                to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            end
            t_temp_landing_m = to_list_m'; clear to_list;
            to_list_k = [];
            for j=1:length(Kx_flights)
                idx = Kx_flights(j);
                to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            end
            t_temp_landing_k = to_list_k'; clear to_list;
        
            % Add spices
            g_id_m = ones(1,length(t_temp_takeoff_m));
            clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);
            g_id_k = ones(1,length(t_temp_takeoff_k));
            clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
            interval_ = [-3 5];
        
            [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v5_cortex(s',t_temp_takeoff_m',t_temp_takeoff_k',t_temp_landing_m',t_temp_landing_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("cortexBatUnit ",num2str(unit)," : takeoff of cortexFlights To M"),strcat("cortexBatUnit ",num2str(unit)," : takeoff of cortexFlights To K"),1,1,Mx_flights,Kx_flights,cortex_flight_struct_resort);
            [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v6_cortex(s',t_temp_takeoff_m',t_temp_takeoff_k',t_temp_landing_m',t_temp_landing_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("cortexBatUnit ",num2str(unit)," : takeoff of cortexFlights To M"),strcat("cortexBatUnit ",num2str(unit)," : takeoff of cortexFlights To K"),1,1,Mx_flights,Kx_flights,cortex_flight_struct_resort,clus);
        else
            Unit_Summary.MvK_sig_pre = 0;
        end

        if mean(sig_land) > 0.89
            s = B_ephys_data.TT_unit(unit).AlignedTimestamps/1e6;
    
            % Take all usec times of flight takeoff CHANGE THIS FOR THE
            % CATEGORIES!!!!
            to_list_m = [];
            for j=1:length(Mx_flights)
                idx = Mx_flights(j);
                to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend_trimmed/120];
            end
            t_temp_takeoff_m = to_list_m'; clear to_list;
            to_list_k = [];
            for j=1:length(Kx_flights)
                idx = Kx_flights(j);
                to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend_trimmed/120];
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
        
            [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v6_cortex_LAND(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("ciholasBatUnit ",num2str(unit)," : landing of ciholasFlights To M"),strcat("ciholasBatUnit ",num2str(unit)," : takeoff of ciholasFlights To K"),1,1,Mx_flights,Kx_flights,cortex_flight_struct_resort,clus);
        else
            Unit_Summary.MvK_sig_land = 0;
        end
    end

    % CONTROLS:
    % Plot the location of the other bat at the time of takeoff
    figure('name','Flights to M (K = yellow), (M = cyan), (OB = green)'); hold on; 
    for j=1:length(Mx_flights)
        i = Mx_flights(j);
        %plot3(cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,1),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,2),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,3),'Color',[0.8 0 0]);
        subplot(2,length(Mx_flights),j); hold on;  title(strcat("Positions at flight ",num2str(Mx_flights(j))," takeoff")); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        plot3(cortex_flight_struct_resort{i}.pos(1,:), cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:),'Color',[0 0 1]);
        
        % Plot human positions
        scatter3(cortex_flight_struct_resort{i}.pos_K(1:10,1), cortex_flight_struct_resort{i}.pos_K(1:10,2),cortex_flight_struct_resort{i}.pos_K(1:10,3),'y','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_M(1:10,1), cortex_flight_struct_resort{i}.pos_M(1:10,2),cortex_flight_struct_resort{i}.pos_M(1:10,3),'cyan','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_ciholas_bat(1:10,1),cortex_flight_struct_resort{i}.pos_ciholas_bat(1:10,2),cortex_flight_struct_resort{i}.pos_ciholas_bat(1:10,3),'green','filled');
        scatter3(cortex_flight_struct_resort{i}.pos(1,1),cortex_flight_struct_resort{i}.pos(2,1),cortex_flight_struct_resort{i}.pos(3,1),'black','filled');
        hold off;
    end
    for j=1:length(Mx_flights)
        i = Mx_flights(j);
        %plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,1),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,2),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,3),'Color',[0.8 0 0]);
        subplot(2,length(Mx_flights),j+length(Mx_flights)); hold on;  title(strcat("Positions at flight ",num2str(Mx_flights(j))," landing")); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        plot3(cortex_flight_struct_resort{i}.pos(1,:), cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:),'Color',[0 0 1]);
        
        % Plot human positions
        scatter3(cortex_flight_struct_resort{i}.pos_K(end-10:end,1), cortex_flight_struct_resort{i}.pos_K(end-10:end,2),cortex_flight_struct_resort{i}.pos_K(end-10:end,3),'y','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_M(end-10:end,1), cortex_flight_struct_resort{i}.pos_M(end-10:end,2),cortex_flight_struct_resort{i}.pos_M(end-10:end,3),'cyan','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_ciholas_bat(end-10:end,1),cortex_flight_struct_resort{i}.pos_ciholas_bat(end-10:end,2),cortex_flight_struct_resort{i}.pos_ciholas_bat(end-10:end,3),'green','filled');
        scatter3(cortex_flight_struct_resort{i}.pos(1,end),cortex_flight_struct_resort{i}.pos(2,end),cortex_flight_struct_resort{i}.pos(3,end),'black','filled');
        hold off;
    end
    hold off;

    % Plot the location of the other bat at the time of takeoff
    figure('name','Flights to K (K = yellow), (M = cyan), (OB = green)'); hold on; 
    for j=1:length(Kx_flights)
        i = Kx_flights(j);
        %plot3(cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,1),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,2),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,3),'Color',[0.8 0 0]);
        subplot(2,length(Kx_flights),j); hold on;  title(strcat("Positions at flight ",num2str(Kx_flights(j))," takeoff")); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        plot3(cortex_flight_struct_resort{i}.pos(1,:), cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:),'Color',[0 0 1]);
        
        % Plot human positions
        scatter3(cortex_flight_struct_resort{i}.pos_K(1:10,1), cortex_flight_struct_resort{i}.pos_K(1:10,2),cortex_flight_struct_resort{i}.pos_K(1:10,3),'y','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_M(1:10,1), cortex_flight_struct_resort{i}.pos_M(1:10,2),cortex_flight_struct_resort{i}.pos_M(1:10,3),'cyan','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_ciholas_bat(1:10,1),cortex_flight_struct_resort{i}.pos_ciholas_bat(1:10,2),cortex_flight_struct_resort{i}.pos_ciholas_bat(1:10,3),'green','filled');
        scatter3(cortex_flight_struct_resort{i}.pos(1,1),cortex_flight_struct_resort{i}.pos(2,1),cortex_flight_struct_resort{i}.pos(3,1),'black','filled');
        hold off;
    end
    for j=1:length(Kx_flights)
        i = Kx_flights(j);
        %plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,1),...
        %                   ciholas_r(ciholas_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,2),...
        %                   cortex_r(cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend,3),'Color',[0.8 0 0]);
        subplot(2,length(Kx_flights),j+length(Kx_flights)); hold on;  title(strcat("Positions at K flight ",num2str(Kx_flights(j))," landing")); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        plot3(cortex_flight_struct_resort{i}.pos(1,:), cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:),'Color',[0 0 1]);
        
        % Plot human positions
        scatter3(cortex_flight_struct_resort{i}.pos_K(end-10:end,1), cortex_flight_struct_resort{i}.pos_K(end-10:end,2),cortex_flight_struct_resort{i}.pos_K(end-10:end,3),'y','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_M(end-10:end,1), cortex_flight_struct_resort{i}.pos_M(end-10:end,2),cortex_flight_struct_resort{i}.pos_M(end-10:end,3),'cyan','filled');
        scatter3(cortex_flight_struct_resort{i}.pos_ciholas_bat(end-10:end,1),cortex_flight_struct_resort{i}.pos_ciholas_bat(end-10:end,2),cortex_flight_struct_resort{i}.pos_ciholas_bat(end-10:end,3),'green','filled');
        scatter3(cortex_flight_struct_resort{i}.pos(1,end),cortex_flight_struct_resort{i}.pos(2,end),cortex_flight_struct_resort{i}.pos(3,end),'black','filled');
        hold off;
    end
    hold off;
end
end