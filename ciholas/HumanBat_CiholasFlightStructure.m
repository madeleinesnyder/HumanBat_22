function [ciholas_flight_struct] = HumanBat_CiholasFlightStructure(exp_data_path,ciholas_flights,cortex_flights,B_ephys_data,redo_with_trimmed_flag)

    %% For each CIHOLAS flight, make structure with all data streams
    %% OR re-do the ephys if flag is up
        % Fields:
        %    fstart: 61227
        %       fend: 61826
        %        pos: [600×3 double]
        %      pos_K: [600×3 double]
        %      pos_M: [600×3 double]
     %pos_cortex_bat: [3×600 double]
        %      fclus: 1
        %      ephys: {1×12 cell}
        %    ephys_t: 4991667
     %fstart_trimmed: 61227
       %fend_trimmed: 61842
      %ephys_trimmed: {1×12 cell}
    %ephys_trimmed_t: 5125001
     %tripod_landing: 8
     %tripod_takeoff: 10

    if redo_with_trimmed_flag == 0
    % Make a structure containg all relevent information for each ciholas flight

    % Load in human and cihiolas bat data 
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));

    cortex_fs = 120;
    for flight_num=1:size(ciholas_flights.flights,1)

         % Get start and stop indexes in cortex data stream
         f_start = ciholas_flights.flights{flight_num,'smp1'};      ciholas_flight_struct{flight_num}.fstart = f_start;
         f_end = ciholas_flights.flights{flight_num,'smp2'};        ciholas_flight_struct{flight_num}.fend = f_end;

         % Get 3d position for that start and stop in ciholas bat data stream
         flightB = ciholas_r(f_start:f_end,:);
         %nanx=isnan(flightB); t=1:numel(flightB); flightB(nanx)=interp1(t(~nanx),flightB(~nanx),t(nanx));   
         %flightB_sm = [smoothdata(flightB(:,1),'gaussian',50),smoothdata(flightB(:,2),'gaussian',50),smoothdata(flightB(:,3),'gaussian',50)];
     
         % Uncomment to see smoothing of flight
         figure(); hold on; plot3(flightB(:,1),flightB(:,2),flightB(:,3));      ciholas_flight_struct{flight_num}.pos = flightB;
     
         % Get 3d position for that start and stop in ciholas human data stream
         cortex_sec_start = f_start/cortex_fs;      ciholas_idx_start = find(round(human_t,4)==round(cortex_sec_start,4)); % Find where in ciholas datastream the seconds match
         cortex_sec_end = f_end/cortex_fs;      ciholas_idx_end = find(round(human_t,4)==round(cortex_sec_end,4)); % Find where in ciholas datastream the seconds match
         flightK = human_r(ciholas_idx_start:ciholas_idx_end,:,3);          ciholas_flight_struct{flight_num}.pos_K = flightK;
         flightM = human_r(ciholas_idx_start:ciholas_idx_end,:,4);          ciholas_flight_struct{flight_num}.pos_M = flightM;
     
         % Get 3d position for that start and stop in cortex bat data stream
         cortex_sec_start = f_start/cortex_fs;      cortex_idx_start = cortex_sec_start*cortex_fs; % Find where in ciholas datastream the seconds match
         cortex_sec_end = f_end/cortex_fs;      cortex_idx_end = cortex_sec_end*cortex_fs; % Find where in ciholas datastream the seconds match
         flightCortexB = cortex_flights.trajectoriesContinous(:,cortex_idx_start:cortex_idx_end)*1000;     ciholas_flight_struct{flight_num}.pos_cortex_bat = flightCortexB;
         
         % Add cluster id from fclus to the struct
         ciholas_flight_struct{flight_num}.fclus = ciholas_flights.flights{flight_num,'fclus'};

         if nansum(nansum(flightCortexB)) ~= 0
            nc=0;
            nanx=isnan(flightCortexB); t=1:numel(flightCortexB); flightCortexB(nanx)=interp1(t(~nanx),flightCortexB(~nanx),t(nanx));   
            flightCortexB_sm = [smoothdata(flightCortexB(:,1),'gaussian',50),smoothdata(flightCortexB(:,2),'gaussian',50),smoothdata(flightCortexB(:,3),'gaussian',50)];    %cortex_flight_struct{flight_num}.pos = flightB_sm;    
         else
             nc = 1;
             disp("No Cortex Data available for this flight");
         end
     
         % Plot the cortex flight, human flights, and ciholas bat flight during that time 
         figure(); hold on; plot3(flightB(:,1),flightB(:,2),flightB(:,3),'r');  
         plot3(flightK(:,1),flightK(:,2),flightK(:,3)); plot3(flightM(:,1),flightM(:,2),flightM(:,3)); 
         if nc==0
             plot3(flightCortexB(1,:),flightCortexB(2,:),flightCortexB(3,:)); 
         end
         title(strcat("Flight ",num2str(flight_num))); hold off;
     
         % flightB - or ciholas_r is the CIHOLAS bat
         % flightCortexB - or cortex_flights.trajectoriesContinuous or Raw is CORTEX
         
        %% Find segment of ephys data corrosponding to that flight 
        % (Right now ephys is only from logger13, which is CIHOLAS/Blondie
        % (Find ephys points closest to the start of the flight according to ciholas data)
        flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
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
                ciholas_flight_struct{flight_num}.ephys{i} = 0;
            else
                disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during flight ",num2str(flight_num)));
                for j=1:length(ephys_idxs_in_range)
                    temp_ts = B_ephys_data.TT_unit(i).AlignedTimestamps(ephys_idxs_in_range(j));
                    temp_rr(round(temp_ts-ephys_f_start)) = 1;
                end
                %rr{flight_num,i} = find(temp_rr==1);
                ciholas_flight_struct{flight_num}.ephys{i} = find(temp_rr==1);
            end
        end

        %Add ephys segment for this flight from each unit to the cortex data struct
        ciholas_flight_struct{flight_num}.ephys_t = length(tt{flight_num});
     
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
        close all;
    end

    % Plot the number of Ciholas spikes during each Ciholas flight for each neuron in box &
    % whisker
    flight_sc=[];
    for j=1:length(B_ephys_data.TT_unit)
        if j==11
            continue
        elseif j==14
            continue
        else
            for i=1:size(ciholas_flights.flights(:,'t1'),1)
                flight_sc(i,j) = length(ciholas_flight_struct{i}.ephys{j});
            end
        end
    end
    % Flight_sc = every row is a flight, every column is a puytative unit
    figure(); hold on; boxplot(flight_sc); xlabel("Unit"); ylabel("Mean Spike Count (101 flights)"); title("Average spike count across all 101 flights per putative unit"); hold off;
        
    % Plot normalizing for time
    for j=1:length(B_ephys_data.TT_unit)
        if j==11
            continue
        elseif j==14
            continue
        else
            for i=1:size(ciholas_flights.flights(:,'t1'),1)
                flight_sc_norm(i,j) = length(ciholas_flight_struct{i}.ephys{j})/((ciholas_flight_struct{i}.fend-ciholas_flight_struct{i}.fstart)/120);
            end
        end
    end
    figure(); hold on; boxplot(flight_sc_norm); xlabel("Unit"); ylabel("Normalized Mean Spike Count (101 flights)"); title("Average spike count across all 101 flights per putative unit"); hold off;
    
    %% Angelo's raster code
    % Make a raster plot for each flight all putative units
    % flight
%     for i=1:12
%         s = B_ephys_data.TT_unit(i).AlignedTimestamps/1e6;
%         
%         % Take all usec times of flight takeoff
%         t_temp_takeoff = ciholas_flights.flights.t1;
%         
%         % Add spices
%         g_id = ones(1,length(t_temp_takeoff));
%         clr = repmat(lines(1),length(t_temp_takeoff),1);
%         interval_ = [-3 5];
% 
%         [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s',t_temp_takeoff',g_id,clr,interval_,strcat("Spikes from Ciholas Bat Unit ",num2str(i)," aligned to takeoff of Ciholas flight"),1,1);
%     end

    %% If we are simply re-doing the ephys, load the pruned struct that
    % already exists
    elseif redo_with_trimmed_flag == 1
        logger=13;
        load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
        load(strcat(exp_data_path,'ciholas/pruned_resorted_ciholas_bat_final_flight_structure.mat'));

        % If already a trimmed field, kill it
        if isfield(ciholas_flight_struct_resort{1},'ephys_trimmed')
            disp("Clearing ephys trimmed data to replace")
            for jj=1:length(ciholas_flight_struct_resort)
                ciholas_flight_struct_resort{jj}.ephys_trimmed = [];
                ciholas_flight_struct_resort{jj}.ephys_trimmed_t = [];   
            end
        end

        % For every flight, redo the ephys entry
        for ob=1:length(ciholas_flight_struct_resort)
            %if ~isfield(ciholas_flight_struct_resort{ob},'fstart_trimmed')
            %    ciholas_flight_struct_resort{ob}.fstart_trimmed = ciholas_flight_struct_resort{ob}.fstart;
            %    ciholas_flight_struct_resort{ob}.fend_trimmed = ciholas_flight_struct_resort{ob}.fend;
            %end

            f_start = ciholas_flight_struct_resort{ob}.fstart_trimmed;
            f_end = ciholas_flight_struct_resort{ob}.fend_trimmed;

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
                    ciholas_flight_struct_resort{ob}.ephys_trimmed{i} = 0;
                else
                    disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during flight ",num2str(ob)));
                    for j=1:length(ephys_idxs_in_range)
                        temp_ts = B_ephys_data.TT_unit(i).AlignedTimestamps(ephys_idxs_in_range(j));
                        temp_rr(round(temp_ts-ephys_f_start)) = 1;
                    end
                    %rr{flight_num,i} = find(temp_rr==1);
                    ciholas_flight_struct_resort{ob}.ephys_trimmed{i} = find(temp_rr==1);
                end
            end
    
            %Add ephys segment for this flight from each unit to the cortex data struct
            ciholas_flight_struct_resort{ob}.ephys_trimmed_t = length(tt{ob});
            ciholas_flight_struct_resort{ob}.ephys_t = length(ciholas_flight_struct_resort{ob}.ephys_t);
        end
        save(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct_resort','-v7.3');
    end
end
