function [ciholas_rest_struct] = HumanBat_CiholasRestStructure(batdate,logger,exp_data_path,plot_flag)

    %% For each CIHOLAS location, make structure with all data streams
    %% Based on ciholas_flight_struct_resort that has been pruned, trimmed, resorted, and tripodded
    % Fields:
            % Ciholas (SELF) 3d position
            % Cortex (non-self) 3d position
            % Ephys (SELF)
            % Ephys (non-self)
            % Audio
            % K position data
            % M position data
            % fstart and fstart_trimmed (start sample of the flight)
            % fend and fend_trimmed (end sample of the flight)
        % Measures:
            % Distance to each human
            % Distanct to the other bat

    % Ground truth of tripods
    tripods = [-2.4,0.56,1.1;
           -1.6,2.07,1.1;
           2.07,1.6,0.9;
           2.26,0.1,0.9;
           2.1,-1.6,1;
           0.36,-2.2,0.9;
           -0.511, -0.265, 0.610;
           -2.545, 1.540, 1.704;
           1.864, -2.197, 1.732;
           1.925, 2.097, 1.742];

    tripods = tripods*1000;
    cutoff=600;
    cortex_fs = 120;

    % Load in the bat position data and human position data  
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
    load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));

    % Load original flight seg data
    load(strcat(exp_data_path,'ciholas/ciholas_bat_final_flight_structure.mat'));
    % Load pruned, trimmed, tripodded flight seg data
    load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));

    %% 1. Make a structure of the "rest" trials

    %    rstart: 61227
    %       rend: 61826
    %        pos: [600×3 double]
    %      pos_K: [600×3 double]
    %      pos_M: [600×3 double]
 %pos_cortex_bat: [3×600 double]
    %      fclus_pre: 1
    %      fclus_pos: 1
    %      ephys: {1×12 cell}
    %    ephys_t: 4991667
    %     tripod: 1

    ciholas_rest_struct = {};
    for i=1:length(ciholas_flight_struct)
        ciholas_rest_struct{i}.rend = ciholas_flight_struct{i}.fstart;
        if i==1
            ciholas_rest_struct{i}.rstart = 1;
        else
            ciholas_rest_struct{i}.rstart = ciholas_flight_struct{i-1}.fend;
        end

        % Get self position
        ciholas_rest_struct{i}.pos = ciholas_r(ciholas_rest_struct{i}.rstart:ciholas_rest_struct{i}.rend,:);

        % Get 3d position for that start and stop in ciholas human data stream
        cortex_sec_start = ciholas_rest_struct{i}.rstart/cortex_fs;      ciholas_idx_start = find(round(human_t,4)==round(cortex_sec_start,4)); % Find where in ciholas datastream the seconds match
        cortex_sec_end = ciholas_rest_struct{i}.rend/cortex_fs;      ciholas_idx_end = find(round(human_t,4)==round(cortex_sec_end,4)); % Find where in ciholas datastream the seconds match
        flightK = human_r(ciholas_idx_start:ciholas_idx_end,:,3);          ciholas_rest_struct{i}.pos_K = flightK;
        flightM = human_r(ciholas_idx_start:ciholas_idx_end,:,4);          ciholas_rest_struct{i}.pos_M = flightM;
     
        % Get 3d position for that start and stop in cortex bat data stream
        cortex_sec_start = ciholas_rest_struct{i}.rstart/cortex_fs;      cortex_idx_start = cortex_sec_start*cortex_fs; % Find where in ciholas datastream the seconds match
        cortex_sec_end = ciholas_rest_struct{i}.rend/cortex_fs;      cortex_idx_end = cortex_sec_end*cortex_fs; % Find where in ciholas datastream the seconds match
        flightCortexB = cortex_flights.trajectoriesContinous(:,cortex_idx_start:cortex_idx_end)*1000;     ciholas_rest_struct{i}.pos_cortex_bat = flightCortexB;

        % Get the clusters of the flight leading up to and following this rest
        if i==1
            ciholas_rest_struct{i}.fclus_pre = NaN;
            ciholas_rest_struct{i}.fclus_post = ciholas_flight_struct{i}.fclus;
        else
            ciholas_rest_struct{i}.fclus_pre = ciholas_flight_struct{i-1}.fclus;
            ciholas_rest_struct{i}.fclus_post = ciholas_flight_struct{i}.fclus;
        end

        % Get the location of that rest spot if available. Plot to be sure
        for j=1:length(ciholas_flight_struct_resort)
            if ciholas_rest_struct{i}.rend == ciholas_flight_struct_resort{j}.fstart
                ciholas_rest_struct{i}.tripod = ciholas_flight_struct_resort{i}.tripod_takeoff;
                ciholas_rest_struct{i}.rend = ciholas_flight_struct_resort{j}.fstart_trimmed;
            else
                ciholas_rest_struct{i}.tripod = 0;
            end
        end
        if ciholas_rest_struct{i}.tripod == 0
            dist_from_tripods = pdist2(tripods, mean(ciholas_rest_struct{i}.pos), 'euclidean');
            tripod_takeoff=0;
             for kk=1:length(dist_from_tripods)
                if dist_from_tripods(kk) < cutoff
                    ciholas_rest_struct{i}.tripod = kk;
                    tripod_takeoff = tripod_takeoff+1;
                end
            end
            if tripod_takeoff==0
                flight_takeoff(i)=11; %other
            elseif tripod_takeoff>1
                flight_takeoff(i)=12; % multiple
            end
        end
        
        % Plot bat and human positions
        if plot_flag
        figure(); hold on; title(strcat("Location of bat before flight ",num2str(i))); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        scatter3(ciholas_rest_struct{i}.pos(:,1),ciholas_rest_struct{i}.pos(:,2),ciholas_rest_struct{i}.pos(:,3),10,'black','filled');
        scatter3(ciholas_rest_struct{i}.pos_K(:,1),ciholas_rest_struct{i}.pos_K(:,2),ciholas_rest_struct{i}.pos_K(:,3),10,'r','filled');
        scatter3(ciholas_rest_struct{i}.pos_M(:,1),ciholas_rest_struct{i}.pos_M(:,2),ciholas_rest_struct{i}.pos_M(:,3),10,'b','filled'); hold off;

        % Comet it
        figure(); hold on; title(strcat("Location of bat before flight ",num2str(i))); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        comet3(ciholas_rest_struct{i}.pos(:,1),ciholas_rest_struct{i}.pos(:,2),ciholas_rest_struct{i}.pos(:,3));
        comet3(ciholas_rest_struct{i}.pos_K(:,1),ciholas_rest_struct{i}.pos_K(:,2),ciholas_rest_struct{i}.pos_K(:,3));
        comet3(ciholas_rest_struct{i}.pos_M(:,1),ciholas_rest_struct{i}.pos_M(:,2),ciholas_rest_struct{i}.pos_M(:,3)); hold off;

        % Plot velocities
        figure(); 
        subplot(3,1,1); hold on; title("Bat Velocity"); plot(diff(ciholas_rest_struct{i}.pos)); hold off;
        subplot(3,1,2); hold on; title("Kevin Velocity"); plot(diff(ciholas_rest_struct{i}.pos_K)); hold off;
        subplot(3,1,3); hold on; title("Madeleine Velocity"); plot(diff(ciholas_rest_struct{i}.pos_M)); hold off;

        % Plot proximity of Madeleine and Kevin to bat
        figure(); kevin_distance = []; mads_distance = [];
        subplot(3,1,1); hold on; title("Bat Velocity"); plot(diff(ciholas_rest_struct{i}.pos)); hold off;
        subplot(3,1,2); hold on; title("Proximity of Kevin to bat"); 
        for j=1:length(ciholas_rest_struct{i}.pos_K); kevin_distance(j) = sqrt(sum((ciholas_rest_struct{i}.pos_K(j,:) - ciholas_rest_struct{i}.pos(j,:)).^2)); end
        plot(kevin_distance/1000); ylim([0 6]); hold off;
        subplot(3,1,3); hold on; title("Proximity of Madeleine to bat"); 
        for j=1:length(ciholas_rest_struct{i}.pos_M); mads_distance(j) = sqrt(sum((ciholas_rest_struct{i}.pos_M(j,:) - ciholas_rest_struct{i}.pos(j,:)).^2)); end
        plot(mads_distance/1000);  ylim([0 6]); hold off;
        end

        % Get the ephys data associated with this stretch of rest
        flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
        f_start_seconds = ciholas_rest_struct{i}.rstart/cortex_fs;  f_end_seconds = ciholas_rest_struct{i}.rend/cortex_fs;   % Get ephys samples in seconds
        for m=1:length(B_ephys_data.TT_unit)
            % Find the closest ephys sample to f_start_seconds*1e6 
            ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
            ephys_f_end = f_end_seconds*1e6;
    
            % Make a vector tt that is the length of the seconds 1e6 in
            % that flight
            tt{i} = [0:(ephys_f_end-ephys_f_start)];   
            temp_rr = zeros(length(tt{i}),1);

            % Fill in the rest of the timestamps with ephys_f_start
            % subtracted for that ephys stretch
            [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(m).AlignedTimestamps & B_ephys_data.TT_unit(m).AlignedTimestamps < ephys_f_end);
            if isempty(ephys_idxs_in_range)
                disp("This unit has no activity during this flight");
                ciholas_rest_struct{i}.ephys{i} = 0;
            else
                disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during flight ",num2str(i)));
                for pp=1:length(ephys_idxs_in_range)
                    temp_ts = B_ephys_data.TT_unit(m).AlignedTimestamps(ephys_idxs_in_range(pp));
                    temp_rr(round(temp_ts-ephys_f_start)) = 1;
                end
                %rr{flight_num,i} = find(temp_rr==1);
                ciholas_rest_struct{i}.ephys{m} = find(temp_rr==1);
            end
        end

        %A dd ephys segment for this flight from each unit to the cortex data struct
        ciholas_rest_struct{i}.ephys_t = length(tt{i});


    % Naive goal: see if there is firing that aligns to the humans moving when the bats are at rest

    if plot_flag
    for ll=1:length(ciholas_rest_struct{i}.ephys)
        figure(); unit=ll;
        subplot(4,2,1); hold on; title(strcat("Ephys train of unit ",num2str(ll)," Ciholas bat at rest after flight ",num2str(i)));
        ephys_ones = zeros(1,ciholas_rest_struct{i}.ephys_t);
        try
            ephys_ones(ciholas_rest_struct{i}.ephys{unit}) = 1;
        catch
            disp("No spikes for this unit")
            continue
        end
        plot(ephys_ones); hold off; kev_distance = []; mads_distance = [];
        subplot(4,2,3); hold on; title("Madeleine Velocity");
        plot(diff(ciholas_rest_struct{i}.pos_M)); hold off;
        subplot(4,2,5); hold on; title("Proximity of Madeleine to bat"); 
        for j=1:length(ciholas_rest_struct{i}.pos_M); mads_distance(j) = sqrt(sum((ciholas_rest_struct{i}.pos_M(j,:) - ciholas_rest_struct{i}.pos(j,:)).^2)); end
        plot(mads_distance/1000);  ylim([0 6]); hold off;
        subplot(4,2,7); hold on; title("Mads flight with bat ephys plotted on top");
        plot3(ciholas_rest_struct{i}.pos_M(:,1),ciholas_rest_struct{i}.pos_M(:,2),ciholas_rest_struct{i}.pos_M(:,3));
        xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        flight_dur = (ciholas_rest_struct{i}.rend -  ciholas_rest_struct{i}.rstart);
        flight_vec = [1:flight_dur]./120;
        if ciholas_rest_struct{i}.ephys{unit} == 0
            disp("No unit firing on this flight");
        else
            for j=1:length(ciholas_rest_struct{i}.ephys{unit})
                nearest_ephys_point = dsearchn(flight_vec',ciholas_rest_struct{i}.ephys{unit}(j)/1e6);
                scatter3(human_r(ciholas_rest_struct{i}.rstart+nearest_ephys_point,1,4),human_r(ciholas_rest_struct{i}.rstart+nearest_ephys_point,2,4),human_r(ciholas_rest_struct{i}.rstart+nearest_ephys_point,3,4),'b','filled');
            end
        end
        scatter3(ciholas_rest_struct{i}.pos(:,1),ciholas_rest_struct{i}.pos(:,2),ciholas_rest_struct{i}.pos(:,3),10,'black','filled');
        hold off;
       
        subplot(4,2,2); hold on; title(strcat("Ephys train of unit ",num2str(ll)," Ciholas bat at rest after flight ",num2str(i)));
        ephys_ones = zeros(1,ciholas_rest_struct{i}.ephys_t);
        ephys_ones(ciholas_rest_struct{i}.ephys{unit}) = 1;
        plot(ephys_ones); hold off;
        subplot(4,2,4); hold on; title("Kevin Velocity");
        plot(diff(ciholas_rest_struct{i}.pos_K)); hold off;
        subplot(4,2,6); hold on; title("Proximity of Kevin to bat"); 
        for j=1:length(ciholas_rest_struct{i}.pos_K); kev_distance(j) = sqrt(sum((ciholas_rest_struct{i}.pos_K(j,:) - ciholas_rest_struct{i}.pos(j,:)).^2)); end
        plot(kev_distance/1000);  ylim([0 6]); hold off;
        subplot(4,2,8); hold on; title("Kev flight with bat ephys plotted on top");
        plot3(ciholas_rest_struct{i}.pos_K(:,1),ciholas_rest_struct{i}.pos_K(:,2),ciholas_rest_struct{i}.pos_K(:,3));
        xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        flight_dur = (ciholas_rest_struct{i}.rend -  ciholas_rest_struct{i}.rstart);
        flight_vec = [1:flight_dur]./120;
        if ciholas_rest_struct{i}.ephys{unit} == 0
            disp("No unit firing on this flight");
        else
            for j=1:length(ciholas_rest_struct{i}.ephys{unit})
                nearest_ephys_point = dsearchn(flight_vec',ciholas_rest_struct{i}.ephys{unit}(j)/1e6);
                scatter3(human_r(ciholas_rest_struct{i}.rstart+nearest_ephys_point,1,3),human_r(ciholas_rest_struct{i}.rstart+nearest_ephys_point,2,3),human_r(ciholas_rest_struct{i}.rstart+nearest_ephys_point,3,3),'b','filled');
            end
        end
        scatter3(ciholas_rest_struct{i}.pos(:,1),ciholas_rest_struct{i}.pos(:,2),ciholas_rest_struct{i}.pos(:,3),10,'black','filled');
        hold off;
    end
    end
    end

    %% Concatenate all the time where the bat is at rest in a given location
    location_struct = {}; 
    x_locations = [1,2,3,4,5,6,7,8,9,10]; colorjet = jet(length(x_locations));
    for i=1:length(x_locations)
        temp_pos = [];
        for j=1:length(ciholas_rest_struct)
            if ciholas_rest_struct{j}.tripod == x_locations(i)
                temp_pos = [temp_pos;ciholas_rest_struct{j}.pos];
            end
        end
        location_struct{i}.pos = temp_pos;
    end

    figure(); hold on; title("All resting positions"); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    for i=1:length(x_locations)
        scatter3(location_struct{i}.pos(:,1),location_struct{i}.pos(:,2),location_struct{i}.pos(:,3),10,colorjet(i,:),'filled'); 
    end

    %% Add ephys to the locations
%     for i=1:length(x_locations)
%         disp(strcat("Working on location ",num2str(i)));
%         for kk = 1:length(ciholas_rest_struct{1}.ephys)
%             temp_ephys = [];
%             for j=1:length(ciholas_rest_struct)
%                 if ciholas_rest_struct{j}.tripod == x_locations(i)
%                     e_bin = zeros(1,ciholas_rest_struct{j}.ephys_t);
%                     if ~isempty(ciholas_rest_struct{j}.ephys{kk}) 
%                         if ciholas_rest_struct{j}.ephys{kk}(1) ~= 0
%                             e_bin(ciholas_rest_struct{j}.ephys{kk}) = 1;
%                         end
%                     end
%                     temp_ephys = [temp_ephys,e_bin];
%                 end
%             end
%             location_struct{i}.ephys{kk} = find(temp_ephys==1);
%             location_struct{i}.ephys_t = length(temp_ephys);
%         end
%     end
% 
%     figure(); ephys_t_vec = plot(location_struct{1}.ephys{8})


    

      


end
