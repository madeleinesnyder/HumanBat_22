function [sig_pre,sig_land] = HumanBat_MvK_Permutation_Clusters(clus,unit,Flight_Group_Matrix,ciholas_flight_struct_resort,B_ephys_data)

% Permutation test for CLUSTERS

% Test whether the difference in firing rate before takeoff, or at
% landing between K and M is significant.

% Each unit is going through 3 tests

%% For all Clusters.... for Ciholas
for clust = clus
    % Find the flights to M and to K
    Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
    Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
    Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
    Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)
    
    % IF there isn't enough data, dn't test
    if length(Mx_flights) < 3 | length(Kx_flights) < 3
        disp(strcat("Not enough flight data for cluster ",num2str(clus)));
        sig_units_pre_clus{clust} = [];
        sig_units_dur_clus{clust} = [];
        sig_units_post_clus{clust} = [];
    else
    
        %% CLUSTERS: Permutation test to see whether the 1s interval before flight is sig
        h_m = []; p_m = []; h_k = []; p_k = []; sec=1;
        for uu=unit %1:length(ciholas_flight_struct_resort{1}.ephys_trimmed)
            mean_spike_rate_pre = [];
            unit = uu;
            All_clus_flights = [Mx_flights,Kx_flights];
            sec=1; cortex_fs = 120;
            for n=1:500
                % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
                % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
                % 3. Take the mean. Repeat
                CFV = ones(length(All_clus_flights),1);
                random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
                CFV(random_indexes_subset1) = 0;
                random_indexes_subset2 = find(CFV==1);
                for i=1:length(random_indexes_subset1)
                    pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart_trimmed-sec*120;
                    pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart_trimmed;
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
                        numSpikes_s1(i) = 0;
                    else
                        for j=1:length(ephys_idxs_in_range)
                            temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
                            temp_rr(round(temp_ts-ephys_f_start)) = 1;
                        end
                        %rr{flight_num,i} = find(temp_rr==1);
                        numSpikes_s1(i) = length(find(temp_rr==1));
                    end
                end
                % Now for subset2
                for i=1:length(random_indexes_subset2)
                    pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart_trimmed-sec*120;
                    pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart_trimmed;
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
                        numSpikes_s2(i) = 0;
                    else
                        for j=1:length(ephys_idxs_in_range)
                            temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
                            temp_rr(round(temp_ts-ephys_f_start)) = 1;
                        end
                        %rr{flight_num,i} = find(temp_rr==1);
                        numSpikes_s2(i) = length(find(temp_rr==1));
                    end
                end
                mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
                mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
                mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
            end
            
            % Calculate the mean spike rate for those Kx flights
            clear numSpikes;
            for i=1:length(Kx_flights)
                pre_flight_start = ciholas_flight_struct_resort{Kx_flights(i)}.fstart_trimmed-sec*120;
                pre_flight_end = ciholas_flight_struct_resort{Kx_flights(i)}.fstart_trimmed;
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
                pre_flight_start = ciholas_flight_struct_resort{Mx_flights(i)}.fstart_trimmed-sec*120;
                pre_flight_end = ciholas_flight_struct_resort{Mx_flights(i)}.fstart_trimmed;
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
    
            % Calculate difference between m and k spike rates (abs?)
            mean_spike_rate_MK_diff = mean_spike_rate_pre_M-mean_spike_rate_pre_K;
    
            % Plot the distribution of the shuffle, with the real value
            figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100); title("Takeoff Distr. of difference in spike rate between flights to M v.s. K")
            
            % Do the ttest
            [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
            clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
        end
        
        % Find significant pre-flight units x CLUSTER
        sig_units_pre_clus{clus} = find(h_m==1);
        if sig_units_pre_clus{clus} == unit
            sig_pre = 1;
        else
            sig_pre = 0;
        end
        
    
         %% CLUSTERS: Permutation test to see whether the Post-flight 1+sec is sig
%         h_m = []; p_m = []; h_k = []; p_k = []; clear mean_spike_rate_pre_K mean_spike_rate_pre_M mean_spike_rate_diff
%         for uu=unit %1:length(ciholas_flight_struct_resort{1}.ephys_trimmed)
%             mean_spike_rate_pre = [];
%             unit = uu;
%             All_clus_flights = [Mx_flights,Kx_flights];
%             sec=1; cortex_fs = 120;
%             for n=1:500
%                 % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
%                 % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
%                 % 3. Take the mean. Repeat
%                 CFV = ones(length(All_clus_flights),1);
%                 random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
%                 CFV(random_indexes_subset1) = 0;
%                 random_indexes_subset2 = find(CFV==1);
%                 for i=1:length(random_indexes_subset1)
%                     pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart_trimmed-sec*120;
%                     pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart_trimmed;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s1(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 % Now for subset2
%                 for i=1:length(random_indexes_subset2)
%                     pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart_trimmed-sec*120;
%                     pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart_trimmed;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s2(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
%                 mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
%                 mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
%             end
%             
%             % Calculate the mean spike rate for those Kx flights
%             clear numSpikes;
%             for i=1:length(Kx_flights)
%                 pre_flight_start = ciholas_flight_struct_resort{Kx_flights(i)}.fend_trimmed;
%                 pre_flight_end = ciholas_flight_struct_resort{Kx_flights(i)}.fend_trimmed+sec*120;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_post_K = mean(numSpikes);
%             
%             % Calculate the mean spike rate for those Mx flights
%             clear numSpikes;
%             for i=1:length(Mx_flights)
%                 pre_flight_start = ciholas_flight_struct_resort{Mx_flights(i)}.fend_trimmed;
%                 pre_flight_end = ciholas_flight_struct_resort{Mx_flights(i)}.fend_trimmed+sec*120;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_post_M = mean(numSpikes);
%     
%             mean_spike_rate_MK_diff = mean_spike_rate_post_M-mean_spike_rate_post_K;
%     
%             %figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100);
%             
%             % Do the ttest
%             [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
%             clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
%         end
%         
%         % Find significant pre-flight units x CLUSTER
%         sig_units_post_clus{clus} = find(h_m==1);
%     
    
         %% CLUSTERS: Permutation test to see whether the -100ms around landing is sig
        h_m = []; p_m = []; h_k = []; p_k = []; sec=1;
        for uu=unit %1:length(ciholas_flight_struct_resort{1}.ephys_trimmed)
            mean_spike_rate_pre = [];
            unit = uu;
            All_clus_flights = [Mx_flights,Kx_flights];
            sec=1; cortex_fs = 120;
            for n=1:500
                % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
                % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
                % 3. Take the mean. Repeat
                CFV = ones(length(All_clus_flights),1);
                random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
                CFV(random_indexes_subset1) = 0;
                random_indexes_subset2 = find(CFV==1);
                for i=1:length(random_indexes_subset1)
                    pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fend_trimmed-0.25*120;
                    pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fend_trimmed+60;
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
                         numSpikes_s1(i) = length(find(temp_rr==1));
                    end
                end
                % Now for subset2
                for i=1:length(random_indexes_subset2)
                    pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fend_trimmed-0.25*120;
                    pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fend_trimmed+60;
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
                         numSpikes_s2(i) = length(find(temp_rr==1));
                    end
                end
                mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
                mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
                mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
            end
            
            % Calculate the mean spike rate for those Kx flights
            clear numSpikes;
            for i=1:length(Kx_flights)
                pre_flight_start = ciholas_flight_struct_resort{Kx_flights(i)}.fend_trimmed-0.25*120;
                pre_flight_end = ciholas_flight_struct_resort{Kx_flights(i)}.fend_trimmed+60;
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
                pre_flight_start = ciholas_flight_struct_resort{Mx_flights(i)}.fend_trimmed-0.25*120;
                pre_flight_end = ciholas_flight_struct_resort{Mx_flights(i)}.fend_trimmed+60;
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
    
            mean_spike_rate_MK_diff = mean_spike_rate_pre_M-mean_spike_rate_pre_K;
    
            figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100); title("Landing Distr. of difference in spike rate between flights to M v.s. K")
            
            % Do the ttest
            [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
            clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
        end
        
        % Find significant pre-flight units x CLUSTER
        sig_units_landing_clus{clus} = find(h_m==1);
        if sig_units_landing_clus{clus} == unit
            sig_land = 1;
        else
            sig_land = 0;
        end
        
     
          %% CLUSTERS: Permutation test to see whether the flight interval is sig
%         h_m = []; p_m = []; h_k = []; p_k = []; sec=1;
%         for uu=unit %1:length(ciholas_flight_struct_resort{1}.ephys_trimmed)
%             mean_spike_rate_pre = [];
%             unit = uu;
%             All_clus_flights = [Mx_flights,Kx_flights];
%             sec=3; cortex_fs = 120;
%             for n=1:500
%                 % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
%                 % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
%                 % 3. Take the mean. Repeat
%                 CFV = ones(length(All_clus_flights),1);
%                 random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
%                 CFV(random_indexes_subset1) = 0;
%                 random_indexes_subset2 = find(CFV==1);
%                 for i=1:length(random_indexes_subset1)
%                     pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart_trimmed;
%                     pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fend_trimmed;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s1(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 % Now for subset2
%                 for i=1:length(random_indexes_subset2)
%                     pre_flight_start = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart_trimmed-sec*120;
%                     pre_flight_end = ciholas_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart_trimmed;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s2(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
%                 mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
%                 mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
%             end
%             
%             % Calculate the mean spike rate for those Kx flights
%             clear numSpikes;
%             for i=1:length(Kx_flights)
%                 pre_flight_start = ciholas_flight_struct_resort{Kx_flights(i)}.fstart_trimmed-sec*120;
%                 pre_flight_end = ciholas_flight_struct_resort{Kx_flights(i)}.fstart_trimmed;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_K = mean(numSpikes);
%             
%             % Calculate the mean spike rate for those Mx flights
%             clear numSpikes;
%             for i=1:length(Mx_flights)
%                 pre_flight_start = ciholas_flight_struct_resort{Mx_flights(i)}.fstart_trimmed-sec*120;
%                 pre_flight_end = ciholas_flight_struct_resort{Mx_flights(i)}.fstart_trimmed;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_M = mean(numSpikes);
%     
%             mean_spike_rate_MK_diff = mean_spike_rate_pre_M-mean_spike_rate_pre_K;
%     
%             %figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100);
%             
%             % Do the ttest
%             [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
%             clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
%         end
%         
%         % Find significant pre-flight units x CLUSTER
%         sig_units_dur_clus{clus} = find(h_m==1);
    end
end





%% For all Clusters.... for Cortex
% for clust = 2:10
%     clus = clust;
%     % Find the flights to M and to K
%     Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
%     Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
%     Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
%     Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)
%     
%     % IF there isn't enough data, dn't test
%     if length(Mx_flights) < 3 | length(Kx_flights) < 3
%         disp(strcat("Not enough flight data for cluster ",num2str(clus)));
%         sig_units_pre_clus{clust} = [];
%         sig_units_dur_clus{clust} = [];
%         sig_units_post_clus{clust} = [];
%     else
%     
%         %% CLUSTERS: Permutation test to see whether the 1s interval before flight is sig
%         h_m = []; p_m = []; h_k = []; p_k = []; sec=1;
%         for uu=1:length(cortex_flight_struct_resort{1}.ephys)
%             mean_spike_rate_pre = [];
%             unit = uu;
%             All_clus_flights = [Mx_flights,Kx_flights];
%             sec=1; cortex_fs = 120;
%             for n=1:500
%                 % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
%                 % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
%                 % 3. Take the mean. Repeat
%                 CFV = ones(length(All_clus_flights),1);
%                 random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
%                 CFV(random_indexes_subset1) = 0;
%                 random_indexes_subset2 = find(CFV==1);
%                 for i=1:length(random_indexes_subset1)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart-sec*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s1(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 % Now for subset2
%                 for i=1:length(random_indexes_subset2)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart-sec*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s2(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
%                 mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
%                 mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
%             end
%             
%             % Calculate the mean spike rate for those Kx flights
%             clear numSpikes;
%             for i=1:length(Kx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Kx_flights(i)}.fstart-sec*120;
%                 pre_flight_end = cortex_flight_struct_resort{Kx_flights(i)}.fstart;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_K = mean(numSpikes);
%             
%             % Calculate the mean spike rate for those Mx flights
%             clear numSpikes;
%             for i=1:length(Mx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Mx_flights(i)}.fstart-sec*120;
%                 pre_flight_end = cortex_flight_struct_resort{Mx_flights(i)}.fstart;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_M = mean(numSpikes);
%     
%             mean_spike_rate_MK_diff = mean_spike_rate_pre_M-mean_spike_rate_pre_K;
%     
%             figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100);
%             
%             % Do the ttest
%             [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
%             clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
%         end
%         
%         % Find significant pre-flight units x CLUSTER
%         sig_units_pre_clus{clus} = find(h_m==1);
%         
%     
%         %% CLUSTERS: Permutation test to see whether the Post-flight 1+sec is sig
%         h_m = []; p_m = []; h_k = []; p_k = []; clear mean_spike_rate_pre_K mean_spike_rate_pre_M mean_spike_rate_diff
%         for uu=1:length(cortex_flight_struct_resort{1}.ephys)
%             mean_spike_rate_pre = [];
%             unit = uu;
%             All_clus_flights = [Mx_flights,Kx_flights];
%             sec=1; cortex_fs = 120;
%             for n=1:500
%                 % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
%                 % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
%                 % 3. Take the mean. Repeat
%                 CFV = ones(length(All_clus_flights),1);
%                 random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
%                 CFV(random_indexes_subset1) = 0;
%                 random_indexes_subset2 = find(CFV==1);
%                 for i=1:length(random_indexes_subset1)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart-sec*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s1(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 % Now for subset2
%                 for i=1:length(random_indexes_subset2)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart-sec*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s2(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
%                 mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
%                 mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
%             end
%             
%             % Calculate the mean spike rate for those Kx flights
%             clear numSpikes;
%             for i=1:length(Kx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Kx_flights(i)}.fend;
%                 pre_flight_end = cortex_flight_struct_resort{Kx_flights(i)}.fend+sec*120;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_post_K = mean(numSpikes);
%             
%             % Calculate the mean spike rate for those Mx flights
%             clear numSpikes;
%             for i=1:length(Mx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Mx_flights(i)}.fend;
%                 pre_flight_end = cortex_flight_struct_resort{Mx_flights(i)}.fend+sec*120;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_post_M = mean(numSpikes);
%     
%             mean_spike_rate_MK_diff = mean_spike_rate_post_M-mean_spike_rate_post_K;
%     
%             figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100);
%             
%             % Do the ttest
%             [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
%             clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
%         end
%         
%         % Find significant pre-flight units x CLUSTER
%         sig_units_post_clus{clus} = find(h_m==1);
%     
%     
%          %% CLUSTERS: Permutation test to see whether the -100ms around landing is sig
%         h_m = []; p_m = []; h_k = []; p_k = [];
%         for uu=1:length(cortex_flight_struct_resort{1}.ephys)
%             mean_spike_rate_pre = [];
%             unit = uu;
%             All_clus_flights = [Mx_flights,Kx_flights];
%             sec=1; cortex_fs = 120;
%             for n=1:500
%                 % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
%                 % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
%                 % 3. Take the mean. Repeat
%                 CFV = ones(length(All_clus_flights),1);
%                 random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
%                 CFV(random_indexes_subset1) = 0;
%                 random_indexes_subset2 = find(CFV==1);
%                 for i=1:length(random_indexes_subset1)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fend-0.25*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fend+60;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s1(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 % Now for subset2
%                 for i=1:length(random_indexes_subset2)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart-sec*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s2(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
%                 mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
%                 mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
%             end
%             
%             % Calculate the mean spike rate for those Kx flights
%             clear numSpikes;
%             for i=1:length(Kx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Kx_flights(i)}.fstart-sec*120;
%                 pre_flight_end = cortex_flight_struct_resort{Kx_flights(i)}.fstart;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_K = mean(numSpikes);
%             
%             % Calculate the mean spike rate for those Mx flights
%             clear numSpikes;
%             for i=1:length(Mx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Mx_flights(i)}.fstart-sec*120;
%                 pre_flight_end = cortex_flight_struct_resort{Mx_flights(i)}.fstart;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_M = mean(numSpikes);
%     
%             mean_spike_rate_MK_diff = mean_spike_rate_pre_M-mean_spike_rate_pre_K;
%     
%             figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100);
%             
%             % Do the ttest
%             [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
%             clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
%         end
%         
%         % Find significant pre-flight units x CLUSTER
%         sig_units_landing_clus{clus} = find(h_m==1);
%         
%      
%          %% CLUSTERS: Permutation test to see whether the flight interval is sig
%         h_m = []; p_m = []; h_k = []; p_k = []; sec=1;
%         for uu=1:length(cortex_flight_struct_resort{1}.ephys)
%             mean_spike_rate_pre = [];
%             unit = uu;
%             All_clus_flights = [Mx_flights,Kx_flights];
%             sec=3; cortex_fs = 120;
%             for n=1:500
%                 % 1. Choose randomly length(Mx_flights) from the pool of all that cluster
%                 % 2. Calculate the spike rate from -3 to 0 at that interval for each flight
%                 % 3. Take the mean. Repeat
%                 CFV = ones(length(All_clus_flights),1);
%                 random_indexes_subset1 = randperm(length(All_clus_flights),length(Mx_flights));
%                 CFV(random_indexes_subset1) = 0;
%                 random_indexes_subset2 = find(CFV==1);
%                 for i=1:length(random_indexes_subset1)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fstart;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset1(i))}.fend;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s1(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 % Now for subset2
%                 for i=1:length(random_indexes_subset2)
%                     pre_flight_start = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart-sec*120;
%                     pre_flight_end = cortex_flight_struct_resort{All_clus_flights(random_indexes_subset2(i))}.fstart;
%                     flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                     f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                     % Find the closest ephys sample to f_start_seconds*1e6 
%                     ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                     ephys_f_end = f_end_seconds*1e6;
%             
%                     % Make a vector tt that is the length of the seconds 1e6 in
%                     % that flight
%                     tt = [0:(ephys_f_end-ephys_f_start)];   
%                     temp_rr = zeros(length(tt),1);
%             
%                     % Fill in the rest of the timestamps with ephys_f_start
%                     % subtracted for that ephys stretch
%                     [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                     if isempty(ephys_idxs_in_range)
%                         %disp("This unit has no activity during this flight");
%                         numSpikes(i) = 0;
%                     else
%                         for j=1:length(ephys_idxs_in_range)
%                             temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                             temp_rr(round(temp_ts-ephys_f_start)) = 1;
%                         end
%                         %rr{flight_num,i} = find(temp_rr==1);
%                          numSpikes_s2(i) = length(find(temp_rr==1));
%                     end
%                 end
%                 mean_spike_rate_pre(n,1) = mean(numSpikes_s1);
%                 mean_spike_rate_pre(n,2) = mean(numSpikes_s2);
%                 mean_spike_rate_diff(n) = mean(numSpikes_s1)-mean(numSpikes_s2);
%             end
%             
%             % Calculate the mean spike rate for those Kx flights
%             clear numSpikes;
%             for i=1:length(Kx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Kx_flights(i)}.fstart-sec*120;
%                 pre_flight_end = cortex_flight_struct_resort{Kx_flights(i)}.fstart;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_K = mean(numSpikes);
%             
%             % Calculate the mean spike rate for those Mx flights
%             clear numSpikes;
%             for i=1:length(Mx_flights)
%                 pre_flight_start = cortex_flight_struct_resort{Mx_flights(i)}.fstart-sec*120;
%                 pre_flight_end = cortex_flight_struct_resort{Mx_flights(i)}.fstart;
%                 flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
%                 f_start_seconds = pre_flight_start/cortex_fs;  f_end_seconds = pre_flight_end/cortex_fs;   % Get ephys samples in seconds
%                 % Find the closest ephys sample to f_start_seconds*1e6 
%                 ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
%                 ephys_f_end = f_end_seconds*1e6;
%             
%                 % Make a vector tt that is the length of the seconds 1e6 in
%                 % that flight
%                 tt = [0:(ephys_f_end-ephys_f_start)];   
%                 temp_rr = zeros(length(tt),1);
%             
%                 % Fill in the rest of the timestamps with ephys_f_start
%                 % subtracted for that ephys stretch
%                 [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
%                 if isempty(ephys_idxs_in_range)
%                     %disp("This unit has no activity during this flight");
%                     numSpikes(i) = 0;
%                 else
%                     for j=1:length(ephys_idxs_in_range)
%                         temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
%                         temp_rr(floor(temp_ts-ephys_f_start)) = 1;
%                     end
%                 end
%                 numSpikes(i) = length(find(temp_rr==1));
%             end
%             mean_spike_rate_pre_M = mean(numSpikes);
%     
%             mean_spike_rate_MK_diff = mean_spike_rate_pre_M-mean_spike_rate_pre_K;
%     
%             figure(); hold on; hist(mean_spike_rate_diff); stem(mean_spike_rate_MK_diff,100);
%             
%             % Do the ttest
%             [h_m(uu),p_m(uu)] = ttest2(mean_spike_rate_MK_diff,mean_spike_rate_diff);
%             clear mean_spike_rate_pre mean_spike_rate_pre_M mean_spike_rate_pre_K
%         end
%         
%         % Find significant pre-flight units x CLUSTER
%         sig_units_dur_clus{clus} = find(h_m==1);
%     end
% end
end
    
    