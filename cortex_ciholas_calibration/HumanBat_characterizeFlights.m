function [flights,posK,posB,posM] = HumanBat_characterizeFlights(exp_data_path, varargin)
%HumanBat_characterizeFlights Summary of this function goes here
%   Parameters
%   ----------
%   exp_data_path
%       Path to directory that contains cortex/ and ciholas/ processed data
%       folders
%
%   IGNORE THIS: 'segmented_trajectories' (Name, Value) [Optional]
%       segmented trajectories data from ImBat_SegmentTrajectories. Saves
%       some time by not reprocessing this.
%   
%   Outputs
%   -------
%   flights
%       data structure containing flight information


% -------------------------- Input parser --------------------------
p = inputParser;
addRequired(p,'exp_data_path');
addOptional(p,'segmented_trajectories',0); % Output from ImBat_SegTrajectories
parse(p,exp_data_path,varargin{:});

segmented_trajectories = p.Results.segmented_trajectories;

% ------------------------------------------------------------------

ciholas2cortex = load('ciholas2cortex_scaling_factors.mat').ciholas2cortex;
ciholas_file = dir(fullfile(exp_data_path,'ciholas/*cdp*.mat'));
cortex_file = dir(fullfile(exp_data_path,'cortex/*track*.mat'));
ciholas_data = load(fullfile(ciholas_file.folder, ciholas_file.name));
cortex_data = load(fullfile(cortex_file.folder, cortex_file.name));

disp("Aligning Ciholas & Cortex")
alignedCiholasCortex = HumanBat_alignCiholasCortex(cortex_data,ciholas_data,ciholas2cortex); % align ciholas with cortex
cortex = alignedCiholasCortex.cortex;
ciholas = alignedCiholasCortex.ciholas;

cortex_fs = cortex.VideoFrameRate;
ciholas_fs = ciholas.fs;

disp("Formatting Tracking Marker Data")
% Format the tracking Marker data
[Location, Location_interp] = ImBat_formatTracking(cortex.Markers);
Location_times = [1:length(cortex.AnalogSignals)];

disp("Segmenting flights")
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

    % Calculate closest distance between trajectories
    BK_dist = flightB - flightK;
    BK_dist = sqrt(BK_dist(:,1).^2 + BK_dist(:,2).^2 + BK_dist(:,3).^2);
    [minInd_BK, minDist_BK] = min(BK_dist);

    BM_dist = flightB - flightM;
    BM_dist = sqrt(BM_dist(:,1).^2 + BM_dist(:,2).^2 + BM_dist(:,3).^2);
    [minInd_BM, minDist_BM] = min(BM_dist);

    
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

    close_dist_thresh = 1000; % threshold for 2 points to be considered close
    flyby_thresh = 800;

    flights{flight_num}.to_human = 0;
    flights{flight_num}.from_human = 0;
    flights{flight_num}.to_kq = 0;
    flights{flight_num}.to_ms = 0;
    flights{flight_num}.from_kq = 0;
    flights{flight_num}.from_ms = 0;
    flights{flight_num}.flyby_kq = 0;
    flights{flight_num}.flyby_ms = 0;
    flights{flight_num}.repr = '';
    if(median(BK_dist(1:3),'omitnan') < close_dist_thresh)
        % flight starts from Kevin
        flights{flight_num}.from_kq = 1;
        flights{flight_num}.from_human = 1;
        flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'fr:kq,');
    elseif(median(BM_dist(1:3),'omitnan') < close_dist_thresh)
        % flight starts from Madeleine
        flights{flight_num}.from_ms = 1;
        flights{flight_num}.from_human = 1;
        flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'fr:ms,');
    end

    if(median(BK_dist(end-2:end),'omitnan') < close_dist_thresh)
        % flight ends at Kevin
        flights{flight_num}.to_kq = 1;
        flights{flight_num}.to_human = 1;
        flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'to:kq,');
    elseif(median(BM_dist(end-2:end),'omitnan') < close_dist_thresh)
        % flight ends at Madeleine
        flights{flight_num}.to_ms = 1;
        flights{flight_num}.to_human = 1;
        flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'to:ms,');
    end

    if(min(BM_dist) < flyby_thresh && ~flights{flight_num}.from_ms && ~flights{flight_num}.to_ms)
        flights{flight_num}.flyby_ms = 1;
        flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'by:ms,');
    end

    if(min(BK_dist) < flyby_thresh && ~flights{flight_num}.from_kq && ~flights{flight_num}.to_kq)
        flights{flight_num}.flyby_kq = 1;
        flights{flight_num}.repr = strcat(flights{flight_num}.repr, 'by:kq,');
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
end

