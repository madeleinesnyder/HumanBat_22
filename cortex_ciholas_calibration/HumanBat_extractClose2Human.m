%function [] = HumanBat_extractClose2Human(DATE)
    % Script to extract epochs / trials where bat is close to human
    
    % load data
    ciholas2cortex = load('ciholas2cortex_scaling_factors.mat').ciholas2cortex;
    ciholas_data = load(strcat('/home/batlab/Desktop/HumanBat/data/14592/processed/',DATE,'/b149f/ciholas/extracted_',DATE,'_cdp_1.mat'));
    cortex_data = load(strcat('/home/batlab/Desktop/HumanBat/data/14592/processed/',DATE,'/b149f/cortex/',DATE,'_14592_tracking_1_track.mat'));
    
    alignedCiholasCortex = HumanBat_alignCiholasCortex(cortex_data,ciholas_data,ciholas2cortex); % align ciholas with cortex
    cortex = alignedCiholasCortex.cortex;
    ciholas = alignedCiholasCortex.ciholas;
    
    cortex_fs = cortex.VideoFrameRate;
    ciholas_fs = ciholas.fs;
    
    
    % Format the tracking Marker data
    [Location, Location_interp] = ImBat_formatTracking(cortex.Markers);
    Location_times = [1:length(cortex.AnalogSignals)];
    
    % Segment the flights
    [segmented_trajectories] = ImBat_SegTrajectories(Location,Location_times);
    bat_pos = segmented_trajectories.trajectories_continuous.';
    
    posK = alignedCiholasCortex.ciholas.kq.pos;
    
    posM = alignedCiholasCortex.ciholas.ms.pos;
    posB = alignedCiholasCortex.cortex.avgMarkerPos;
    posB = HumanBat_interpolate_nans(posB);
    
    % Bat speed
    speed_x = diff(posB(:,1))*120;
    speed_y = diff(posB(:,2))*120;
    speed_z = diff(posB(:,3))*120;
    speed = sqrt(speed_x.^2 + speed_y.^2 + speed_z.^2);
    speed = medfilt1(speed,5);
    
    % Plot KQMS positions
    figure;
    scatter(posK(:,1), posK(:,2),5,'red','filled');
    hold on;
    scatter(posM(:,1), posM(:,2),5,'blue','filled');
    legend('KQ', 'MS');
    title('KQ & MS positions')
    
    % Plot Bat&KQMS positions
    figure;
    scatter(posK(:,1), posK(:,2),5,'red','filled');
    hold on;
    scatter(posM(:,1), posM(:,2),5,'blue','filled');
    scatter(posB(:,1), posB(:,2),5,'black','filled');
    legend('KQ', 'MS','Bat');
    title('Bat & KQ & MS positions')
    
    
    % For each flight, get segment close to human
    figure;
    tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'tight');
    flights = {}; % Data structure holding information about each flight{i}
    for flight_num = 1:length(segmented_trajectories.flight_starts_times)
        f_start = segmented_trajectories.flight_starts_indx(flight_num);
        f_end = segmented_trajectories.flight_ends_indx(flight_num);
        
        flightK = posK(f_start:f_end,:);
        flightM = posM(f_start:f_end,:);
        flightB = posB(f_start:f_end,:);
        
        % Get start and end positions of flight. Median of 3 data points for
        % more robust estimation.
        startPosB = median(flightB(1:3,:), 'omitnan');
        endPosB = median(flightB(end-2:end,:), 'omitnan');
        
        % Start and end position of humans during flight (most likely will be
        % stationary)
        startPosK = median(flightK(1:3,:),'omitnan');
        endPosK = median(flightK(end-2:end,:),'omitnan');
        startPosM = median(flightM(1:3,:),'omitnan');
        endPosM = median(flightM(end-2:end,:),'omitnan');
    
        close_dist_thresh= 750; % threshold for 2 points to be considered close
     
        flights{flight_num}.to_human = 0;
        flights{flight_num}.from_human = 0;
        flights{flight_num}.to_kq = 0;
        flights{flight_num}.to_ms = 0;
        flights{flight_num}.from_kq = 0;
        flights{flight_num}.from_ms = 0;
        flights{flight_num}.repr = '';
        if(norm(startPosB - startPosK) < close_dist_thresh)
            % flight starts from Kevin
            flights{flight_num}.from_kq = 1;
            flights{flight_num}.from_human = 1;
            flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'from kq, ');
        elseif(norm(startPosB - startPosM) < close_dist_thresh)
            % flight starts from Madeleine
            flights{flight_num}.from_ms = 1;
            flights{flight_num}.from_human = 1;
            flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'from ms, ');
        end
    
        if(norm(endPosB - endPosK) < close_dist_thresh)
            % flight ends at Kevin
            flights{flight_num}.to_kq = 1;
            flights{flight_num}.to_human = 1;
            flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'to kq, ');
        elseif(norm(endPosB - endPosM) < close_dist_thresh)
            % flight ends at Madeleine
            flights{flight_num}.to_ms = 1;
            flights{flight_num}.to_human = 1;
            flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'to ms, ');
        end
    
        
        % Plot Human Bat positions during bat flights
        nexttile;
        scatter(flightK(:,1), flightK(:,2),5,'black','filled');
        hold on;
        scatter(flightM(:,1), flightM(:,2),5,'green','filled');
        scatter(flightB(:,1), flightB(:,2),5,speed(f_start:f_end),'filled');
        title(flights{flight_num}.repr);
        xlim([-3000 3000])
        ylim([-3000 3000])
    end
    nexttile;
    legend('KQ', 'MS','Bat');
    sgtitle('Bat & KQ & MS positions during flights')

%end
