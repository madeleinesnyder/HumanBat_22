function [cortex_flight_struct] = HumanBat_CortexFlightStructure(exp_data_path,cortex_flights,B_ephys_data,redo_with_trimmed_flag)

    %% For each CORTEX flight, make structure with all data streams
    % Fields:
        % Cortex (SELF) 3d position
        % cortex (non-self) 3d position
        % Ephys (SELF)
        % Ephys (non-self)
        % Audio
        % K position data
        % M position data
    % Measures:
        % Distance to each human
        % Distanct to the other bat

    % Make a structure containg all relevent information for each cortex flight

    if redo_with_trimmed_flag == 0

    % Load in human and cihiolas bat data 
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));

    cortex_fs = 120;
    cortex_flight_struct = {};
    for flight_num=1:length(cortex_flights.flight_starts_idx)
        disp(strcat("Processing Flight ",num2str(flight_num)));

        % Get start and stop indexes in cortex data stream
        f_start = cortex_flights.flight_starts_idx(flight_num);        cortex_flight_struct{flight_num}.fstart = f_start;
        f_end = cortex_flights.flight_ends_idx(flight_num);            cortex_flight_struct{flight_num}.fend = f_end;
        
        % Get 3d position for that start and stop in cortex bat data stream
        flightB = cortex_flights.trajectoriesContinous(:,f_start:f_end)*1000';     cortex_flight_struct{flight_num}.pos = flightB;
        nanx=isnan(flightB); t=1:numel(flightB); flightB(nanx)=interp1(t(~nanx),flightB(~nanx),t(nanx));   
        flightB_sm = [smoothdata(flightB(1,:),'gaussian',50);smoothdata(flightB(2,:),'gaussian',50);smoothdata(flightB(3,:),'gaussian',50)];    cortex_flight_struct{flight_num}.pos = flightB_sm;    
    
        % Uncomment to see smoothing of flight
        %figure(); hold on; plot3(flightB_sm(:,1),flightB_sm(:,2),flightB_sm(:,3));   plot3(flightB(:,1),flightB(:,2),flightB(:,3));
    
        % Get 3d position for that start and stop in ciholas human data stream
        cortex_sec_start = f_start/cortex_fs;      ciholas_idx_start = find(round(human_t,4)==round(cortex_sec_start,4)); % Find where in ciholas datastream the seconds match
        cortex_sec_end = f_end/cortex_fs;      ciholas_idx_end = find(round(human_t,4)==round(cortex_sec_end,4)); % Find where in ciholas datastream the seconds match
        flightK = human_r(ciholas_idx_start:ciholas_idx_end,:,3);       cortex_flight_struct{flight_num}.pos_K = flightK;
        flightM = human_r(ciholas_idx_start:ciholas_idx_end,:,4);       cortex_flight_struct{flight_num}.pos_M = flightM;
    
        % Get 3d position for that start and stop in ciholas BAT data stream
        cortex_sec_start = f_start/cortex_fs;      diff_cc = abs(cortex_sec_start-ciholas_t); ciholas_idx_start = find(diff_cc==min(diff_cc)); %ciholas_idx_start = find(ciholas_t==cortex_sec_start); % Find where in ciholas datastream the seconds match
        cortex_sec_end = f_end/cortex_fs;      diff_cc = abs(cortex_sec_end-ciholas_t); ciholas_idx_end = find(diff_cc==min(diff_cc)); %ciholas_idx_end = find(ciholas_t==cortex_sec_end); % Find where in ciholas datastream the seconds match
        flightCiholasB = ciholas_r(ciholas_idx_start:ciholas_idx_end,:);    cortex_flight_struct{flight_num}.pos_ciholas_bat = flightCiholasB;
    
        % Plot the cortex flight, human flights, and ciholas bat flight during that time 
%         figure(); hold on; plot3(flightB_sm(1,:),flightB_sm(2,:),flightB_sm(3,:),'r');   plot3(flightB(1,:),flightB(2,:),flightB(3,:));
%         plot3(flightK(:,1),flightK(:,2),flightK(:,3)); plot3(flightM(:,1),flightM(:,2),flightM(:,3)); 
%         plot3(flightCiholasB(:,1),flightCiholasB(:,2),flightCiholasB(:,3)); 
%         title(strcat("Flight ",num2str(flight_num))); hold off;
    
        % flightCiholasB - or ciholas_r is the CIHOLAS bat
        % flightB - or cortex_flights.trajectoriesContinuous or Raw is CORTEX
        
        %% Find segment of ephys data corrosponding to that flight 
        % (Right now ephys is only from logger13, which is CIHOLAS/Blondie
        % (Find ephys points closest to the start of the flight according to ciholas data)
        flight_ephys_B = {}; ephys=0;   tt={}; rr={};
        f_start_seconds = f_start/cortex_fs;  f_end_seconds = f_end/cortex_fs;   % Get ephys samples in seconds
        for i=1:length(B_ephys_data.TT_unit)
            % Find the closest ephys sample to f_start_seconds*1e6 
            ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
            ephys_f_end = f_end_seconds*1e6;
    
            % Make a vector tt that is the length of the seconds 1e6 in
            % that flight
            tt{flight_num} = [0:(ephys_f_end-ephys_f_start)];   
            temp_rr = zeros(length(tt{flight_num}),1);

            % Fill in the rest of the timestamps with ephys_f_start
            % subtracted for that ephys stretch
            [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(i).AlignedTimestamps & B_ephys_data.TT_unit(i).AlignedTimestamps < ephys_f_end);
            if isempty(ephys_idxs_in_range)
                disp("This unit has no activity during this flight");
                cortex_flight_struct{flight_num}.ephys{i} = 0;
            else
                disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during flight ",num2str(flight_num)));
                for j=1:length(ephys_idxs_in_range)
                    temp_ts = B_ephys_data.TT_unit(i).AlignedTimestamps(ephys_idxs_in_range(j));
                    temp_rr(round(temp_ts-ephys_f_start)) = 1;
                end
                %rr{flight_num,i} = find(temp_rr==1);
                cortex_flight_struct{flight_num}.ephys{i} = find(temp_rr==1);
            end
        end

        % Add ephys segment for this flight from each unit to the cortex data struct
        %cortex_flight_struct{flight_num}.ephys = rr{flight_num,:};
        cortex_flight_struct{flight_num}.ephys_t = length(tt{flight_num});%tt{flight_num};

        % Assign cluster
        for jj=1:length(cortex_flights.clusterIndex)
            if ismember(flight_num,cortex_flights.clusterIndex{jj})
                cortex_flight_struct{flight_num}.fclus = jj;
            end
        end
         
        % Plot spikes ON the flight for each flight
%         colormap_ = jet(12);
%         % 1. Plot the flight itself
%         figure(); hold on;
%         plot3(flightB_sm(1,:),flightB_sm(2,:),flightB_sm(3,:),'black');
%         for i=1:12
%             % 2. Plot the spikes as little red dots along the flightpath
%             temp_spikes_ts = rr{flight_num,i};
%             num_temp_spikes = length(temp_spikes_ts);
%             for j=1:num_temp_spikes
%                 temp_idx = temp_spikes_ts(j);
%                 temp_idx_s = temp_idx/1e6*120;
%                 temp_idx_s = round(temp_idx_s);     if temp_idx_s(1)==0; temp_idx_s(1)=1; end
%                 scatter3(flightB_sm(1,round(temp_idx_s)),flightB_sm(2,round(temp_idx_s)),flightB_sm(3,round(temp_idx_s)),30,colormap_(i,:),'filled');
%             end
%         end
%         title(strcat("Ciholas Ephys activity for Cortex flight ",num2str(flight_num)));
%         hold off;

        disp("Add Audio Data!");
    end
    clear temp_ts temp_rr

    % Plot the number of Ciholas spikes during each Cortex flight for each neuron in box &
    % whisker
    flight_sc=[];
    for j=1:length(B_ephys_data.TT_unit)
        for i=1:length(cortex_flights.flight_starts_idx)
            flight_sc(i,j) = length(cortex_flight_struct{i}.ephys{j});
        end
    end
    % Flight_sc = every row is a flight, every column is a puytative unit
    figure(); hold on; boxplot(flight_sc); xlabel("Unit"); ylabel("Mean Spike Count (101 flights)"); title("Average spike count across all 101 flights per putative unit"); hold off;
        
    % Plot normalizing for time
    for j=1:length(B_ephys_data.TT_unit)
        for i=1:length(cortex_flights.flight_starts_idx)
            flight_sc_norm(i,j) = length(cortex_flight_struct{i}.ephys{j})/((cortex_flight_struct{i}.fend-cortex_flight_struct{i}.fstart)/120);
        end
    end
    figure(); hold on; boxplot(flight_sc_norm); xlabel("Unit"); ylabel("Normalized Mean Spike Count (101 flights)"); title("Average spike count across all 101 flights per putative unit"); hold off;    

    %% IF REDOING FOR TRIMMING!
    elseif redo_with_trimmed_flag == 1
        logger=13;
        load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
        load(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'));

        % If already a trimmed field, kill it
        if isfield(cortex_flight_struct_resort{1},'ephys_trimmed')
            disp("Clearing ephys trimmed data to replace")
            for jj=1:length(cortex_flight_struct_resort)
                cortex_flight_struct_resort{jj}.ephys_trimmed = [];
                cortex_flight_struct_resort{jj}.ephys_trimmed_t = [];   
            end
        end

        % For every flight, redo the ephys entry
        for ob=1:length(cortex_flight_struct_resort)
            %if ~isfield(cortex_flight_struct_resort{ob},'fstart_trimmed')
            %    ciholas_flight_struct_resort{ob}.fstart_trimmed = ciholas_flight_struct_resort{ob}.fstart;
            %    ciholas_flight_struct_resort{ob}.fend_trimmed = ciholas_flight_struct_resort{ob}.fend;
            %end

            f_start = cortex_flight_struct_resort{ob}.fstart_trimmed;
            f_end = cortex_flight_struct_resort{ob}.fend_trimmed;

            flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
            f_start_seconds = f_start/120;  f_end_seconds = f_end/120;   % Get ephys samples in seconds
            for i=1:length(B_ephys_data.TT_unit)
                % Find the closest ephys sample to f_start_seconds*1e6 
                ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
                ephys_f_end = f_end_seconds*1e6;
        
                % Make a vector tt that is the length of the seconds 1e6 in
                % that flight
                tt{ob} = [0:(ephys_f_end-ephys_f_start)];   
                temp_rr = zeros(length(tt{ob}),1);
    
                % Fill in the rest of the timestamps with ephys_f_start
                % subtracted for that ephys stretch
                [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(i).AlignedTimestamps & B_ephys_data.TT_unit(i).AlignedTimestamps < ephys_f_end);
                if isempty(ephys_idxs_in_range)
                    disp("This unit has no activity during this flight");
                    cortex_flight_struct_resort{ob}.ephys_trimmed{i} = 0;
                else
                    disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during flight ",num2str(ob)));
                    for j=1:length(ephys_idxs_in_range)
                        temp_ts = B_ephys_data.TT_unit(i).AlignedTimestamps(ephys_idxs_in_range(j));
                        temp_rr(round(temp_ts-ephys_f_start)) = 1;
                    end
                    %rr{flight_num,i} = find(temp_rr==1);
                    cortex_flight_struct_resort{ob}.ephys_trimmed{i} = find(temp_rr==1);
                end
            end
    
            %Add ephys segment for this flight from each unit to the cortex data struct
            cortex_flight_struct_resort{ob}.ephys_trimmed_t = length(tt{ob});
            cortex_flight_struct_resort{ob}.ephys_t = length(cortex_flight_struct_resort{ob}.ephys_t);
        end
        save(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');
    end

end

