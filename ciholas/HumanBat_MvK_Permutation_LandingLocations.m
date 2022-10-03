% Permutation tests!

%% LOCATIONS: Permutation test to see whether the -3s interval before flight is sig
% diff for Kevin versus Madeleine
for ff=1:12

    flights_to_smush = strcat('_',num2str(ff));  % i.e. '4_' means "from 4", focusing on takeoff
    
    % If it is a "to" flight...
    Location_Flight_List = [];
    if flights_to_smush(1) == '_'
        for j=1:12
            if str2double(flights_to_smush(2:end)) == j
                index_ = j
                disp(strcat("Looking at rasters of takeoffs of all flights TO location ",num2str(index_)," Madeleine versus Kevin"));
            end
        end 
        % Find all the flights in this Location subset
        for i=1:size(ciholas_flight_struct_resort,2)
            if ciholas_flight_struct_resort{i}.tripod_landing == index_
                Location_Flight_List = [Location_Flight_List,i];
            end
        end
    elseif flights_to_smush(2) == '_'
        for j=1:12
            if str2double(flights_to_smush(1:end-1)) == j
                index_ = j;
                disp(strcat("Looking at rasters of takeoffs of all flights FROM location ",num2str(index_)," Madeleine versus Kevin"));
            end
        end
        % Find all the flights in this Location subset
        for i=1:size(ciholas_flight_struct_resort,2)
            if ciholas_flight_struct_resort{i}.tripod_takeoff == index_
                Location_Flight_List = [Location_Flight_List,i];
            end
        end
    end
    Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Location_Flight_List);
    Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
    Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Location_Flight_List);
    Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)
    
    if length(Mx_flights) < 3 | length(Kx_flights) < 3
        disp(strcat("Not enough flight data for location ",num2str(flights_to_smush)));
        sig_units_pre_loc{ff} = [];
        sig_units_dur_loc{ff} = [];
        sig_units_post_loc{ff} = [];
        continue;
    else
    
        %% CLUSTERS: Permutation test to see whether the +3s interval after flight is sig
        h_m = []; p_m = []; h_k = []; p_k = [];
        for uu=1:length(ciholas_flight_struct_resort{1}.ephys)
            mean_spike_rate_pre = [];
            unit = uu;
            All_clus_flights = [Mx_flights,Kx_flights];
            sec=3; cortex_fs = 120;
            for n=1:500
                % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
                % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
                % 3. Take the mean. Repeat
                random_indexes = randperm(length(All_clus_flights),length(Mx_flights));
                for i=1:length(random_indexes)
                    pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes(i))}.fend_trimmed;
                    pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes(i))}.fend_trimmed+sec*120;
                    flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
                    f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
                    % Find the closest ephys sample to f_start_seconds*1e6 
                    ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
                    ephys_f_end = f_end_seconds*1e6;
            
                    % Make a vector tt that is the length of the seconds 1e6 in
                    % that flight
                    tt = [0:(ephys_f_end-ephys_f_start)];   
                    temp_rr = zeros(length(tt),1);
            
                    % Fill in the rest of the timestamps with ephys_f_start
                    % subtracted for that ephys stretch
                    [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
                    if isempty(ephys_idxs_in_range)
                        %disp("This unit has no activity during this flight");
                        numSpikes(i) = 0;
                    else
                        for j=1:length(ephys_idxs_in_range)
                            temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
                            temp_rr(round(temp_ts-ephys_f_start)) = 1;
                        end
                        %rr{flight_num,i} = find(temp_rr==1);
                         numSpikes(i) = length(find(temp_rr==1));
                    end
                end
                mean_spike_rate_pre(n) = mean(numSpikes);
            end
            
            % Calculate the mean spike rate for those Kx flights
            clear numSpikes;
            for i=1:length(Kx_flights)
                pre_flight_start = ciholas_flight_struct_resort{Kx_flights(i)}.fend_trimmed;
                pre_flight_end = ciholas_flight_struct_resort{Kx_flights(i)}.fend_trimmed+sec*120;
                flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
                f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
                % Find the closest ephys sample to f_start_seconds*1e6 
                ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
                ephys_f_end = f_end_seconds*1e6;
            
                % Make a vector tt that is the length of the seconds 1e6 in
                % that flight
                tt = [0:(ephys_f_end-ephys_f_start)];   
                temp_rr = zeros(length(tt),1);
            
                % Fill in the rest of the timestamps with ephys_f_start
                % subtracted for that ephys stretch
                [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
                if isempty(ephys_idxs_in_range)
                    %disp("This unit has no activity during this flight");
                    numSpikes(i) = 0;
                else
                    for j=1:length(ephys_idxs_in_range)
                        temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
                        temp_rr(floor(temp_ts-ephys_f_start)) = 1;
                    end
                end
                numSpikes(i) = length(find(temp_rr==1));
            end
            mean_spike_rate_pre_K = mean(numSpikes);
            
            % Calculate the mean spike rate for those Mx flights
            clear numSpikes;
            for i=1:length(Mx_flights)
                pre_flight_start = ciholas_flight_struct_resort{Mx_flights(i)}.fend_trimmed;
                pre_flight_end = ciholas_flight_struct_resort{Mx_flights(i)}.fend_trimmed+sec*120;
                flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
                f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
                % Find the closest ephys sample to f_start_seconds*1e6 
                ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
                ephys_f_end = f_end_seconds*1e6;
            
                % Make a vector tt that is the length of the seconds 1e6 in
                % that flight
                tt = [0:(ephys_f_end-ephys_f_start)];   
                temp_rr = zeros(length(tt),1);
            
                % Fill in the rest of the timestamps with ephys_f_start
                % subtracted for that ephys stretch
                [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
                if isempty(ephys_idxs_in_range)
                    %disp("This unit has no activity during this flight");
                    numSpikes(i) = 0;
                else
                    for j=1:length(ephys_idxs_in_range)
                        temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
                        temp_rr(floor(temp_ts-ephys_f_start)) = 1;
                    end
                end
                numSpikes(i) = length(find(temp_rr==1));
            end
            mean_spike_rate_pre_M = mean(numSpikes);
            
            % Do the ttest
            [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_pre_M,mean_spike_rate_pre);
            [h_k(uu),p_k(uu)] = ttest2(mean_spike_rate_pre_K,mean_spike_rate_pre);
            clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
        end
        
        % Find significant pre-flight units x CLUSTER
        sig_units_pre_loc_land{ff} = find(h_m==1 | h_k ==1);
        
        %% Rasters of To Madeleine vs Kevin x Cluster (clus); Takeoff 
        for i=1:length(sig_units_pre_loc_land{ff}) %length(B_ephys_data.TT_unit)
            ii=sig_units_pre_loc_land{ff}(i);%cool_list(i);
            s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
            
            % Take all usec times of flight takeoff CHANGE THIS FOR THE
            % CATEGORIES!!!!
            to_list_m = [];
            for j=1:length(Mx_flights)
                idx = Mx_flights(j);
                to_list_m = [to_list_m,ciholas_flight_struct_resort{idx}.fend_trimmed/120];
            end
            t_temp_takeoff_m = to_list_m'; clear to_list;
            to_list_k = [];
            for j=1:length(Kx_flights)
                idx = Kx_flights(j);
                to_list_k = [to_list_k,ciholas_flight_struct_resort{idx}.fend_trimmed/120];
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
        
            [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s',t_temp_takeoff_m',t_temp_takeoff_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("CiholasBatUnit ",num2str(ii)," : takeoff of CiholasFlights To M"),strcat("CiholasBatUnit ",num2str(ii)," : takeoff of CiholasFlights To K"),1,1);
            h = findobj('Type','figure');
            h(1).Position = [10 10 1000 150];
            figure_name = strcat(num2str(batdate),'_14650_U',num2str(ii),'_CiholasTakeoff_MvK_PostFLight_LandingLoc',num2str(ff));
            saveas(h(1),strcat(exp_data_path(1:48),'Figures/by_neuron/',figure_name,'.svg'));
            close all;
        end
        
        %% Plot Kevin versus Madeleine flights in colors
        figure(); hold on;
        for k=1:length(Kx_flights)
            flightnum = Kx_flights(k);
            plot3(ciholas_flight_struct_resort{flightnum}.pos(:,1),ciholas_flight_struct_resort{flightnum}.pos(:,2),ciholas_flight_struct_resort{flightnum}.pos(:,3),'r');
            scatter3(ciholas_flight_struct_resort{flightnum}.pos(1,1),ciholas_flight_struct_resort{flightnum}.pos(1,2),ciholas_flight_struct_resort{flightnum}.pos(1,3),'green','filled');
        end
        for k=1:length(Mx_flights)
            flightnum = Mx_flights(k);
            plot3(ciholas_flight_struct_resort{flightnum}.pos(:,1),ciholas_flight_struct_resort{flightnum}.pos(:,2),ciholas_flight_struct_resort{flightnum}.pos(:,3),'b');
            scatter3(ciholas_flight_struct_resort{flightnum}.pos(1,1),ciholas_flight_struct_resort{flightnum}.pos(1,2),ciholas_flight_struct_resort{flightnum}.pos(1,3),'green','filled');
        end
        xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300])
        title(strcat("All flights: K (red) Flights to M (blue)"));
        hold off; 
        h = findobj('Type','figure');
        figure_name = strcat(num2str(batdate),'_14650_MvK_LandingLocation_',num2str(ff));
        saveas(h(1),strcat(exp_data_path(1:48),'Figures/by_flight/',figure_name,'.svg'));
        close all;

    end
end

