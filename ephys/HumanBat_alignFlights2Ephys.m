function [flights] = HumanBat_alignFlights2Ephys(exp_data_path, date, varargin)

% This is the effective entrypoint for any HumanBat Analysis
    % 1. Process the cortex, ciholas, and ephys data using the Luigi pipeline 
        % a. conda activate humansniffbat
        % b. Navigate to HumanBat dir with humanbat.py script
        % c. python humanbat.py --batname = '2992814650' --batdate='220411'
        % d. To change what processes are run, open Atom to change humanbat.py
    % 2. Make sure all files are in their proper folders (follow the format for 220407)
    % 3. Align cortex, ciholas, and ephys data
        % a. This script calls those scripts!
        % b. This script runs checks on that alignment!
    % 4. Make ciholas and cortex flight structs
        % a. This script calls those scripts!
    % 5. Prune the flight structs of the human-carrying flights and bad flights
        % a. HumanBat_Ciholas_Flights_Ephys_Analysis.m
    % 6. Resort the flights according to cluster
        % a. HumanBat_Ciholas_Flights_Ephys_Analysis.m
    % 7. Trim the flights start and stop times
        % a. HumanBat_Flight_Trimming
    % 8. Get landing and takeoff tripod locations
        % a. HumanBat_LandingLocationTakeoffLocation.m
    % 9. Re-do the ephys --> ephys_trimmed to account for trimming
        % a. [ciholas_flight_struct] = HumanBat_CiholasFlightStructure(exp_data_path,ciholas2cortex,ciholas2cortex,B_ephys_data,1);
    % 10. If you want to look at Spatial Information and significance:
        %. HumanBat_SI.m and HumanBat_SI_red_dots.m
    % 11. If you want to look at rasters of M or K x Cluster or x Location
        % HumanBat_Ciholas_Flights_Ephys_Analysis.m
    % 12. If you want to look at movies aligned to spikes in auditory space:
        % HumanBat_shell_brainMovies(220407,13,5,"top",490,180);

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

%% TEMPORARY:
batdate=220418;
if batdate==220422
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/14650/processed/',num2str(batdate),'/');
else
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
end
logger=15;

% Paths
% ------------------------------------------------------------------
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ImBat'));
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/cortex_ciholas_calibration'));
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ciholas'));
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/cortex'));

% Load Data
% ------------------------------------------------------------------

disp("Loading Cortex, Ciholas, and Ephys Data");
ciholas2cortex = load('ciholas2cortex_scaling_factors.mat').ciholas2cortex;
ciholas_file = dir(fullfile(exp_data_path,'ciholas/*cdp*.mat'));
cortex_file = dir(fullfile(exp_data_path,'cortex/*track*'));
B_ephys_file = dir(fullfile(strcat(exp_data_path,'/ephys/logger',num2str(logger),'/extracted_data/SingleUnit*')));      
ciholas_data = load(fullfile(ciholas_file.folder, ciholas_file.name));
cortex_data = load(fullfile(cortex_file.folder, cortex_file.name));
B_ephys_data = load(strcat(B_ephys_file(end).folder,'/',B_ephys_file(end).name,'/SingleUnits_',num2str(batdate),'.mat'));   
ephys_TTL = load(strcat(exp_data_path,'/ephys/logger',num2str(logger),'/extracted_data/00000_','20',num2str(batdate),'_EVENTS.mat'));
cortex_file = []; ciholas_file =[]; B_ephys_file=[]; 

% Get sampling rates for all data systems
cortex_fs = cortex_data.AnalogFrameRate;
%ciholas_fs = ciholas_data.CDPmtdata.Fs;
ephys_fs = B_ephys_data.TT_unit(1).Fs;

%% CORTEX
% ------------------------------------------------------------------

% Check if Cortex data has been formatted and saved
if exist(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'))
    disp("Cortex Data has been formatted, segmented, and clustered. Loading now...");
    load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
    load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
else
    disp("Formatting Cortex Tracking Marker Data")
    % Format the tracking Marker data
    [Location, Location_interp] = ImBat_formatTracking(cortex_data.Markers);
    Location_times = [1:length(cortex_data.AnalogSignals)];
    
    disp("Segmenting Cortex Flights")
    % Segment the flights
    [segmented_trajectories,cortex_ttl_seconds] = ImBat_SegTrajectories(Location,Location_times,cortex_data);
    %close all;  
    save(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'),'segmented_trajectories','cortex_ttl_seconds');
    Location = []; Location_interp = []; Location_times = [];
    %%
    disp("Clustering Cortex Flights")
    % Cluster the flights (220406 dist = 1.2, splines = 7)
    if batdate==220331
        dist_ = 2.2;
    elseif batdate==220401
        dist_ = 2;
    elseif batdate==220404
        dist_ = 1.4;
    elseif batdate==220406
        dist_ = 1.6;
    elseif batdate==220407
        dist_ = 1.1;
    elseif batdate==220408
        dist_ = 1.6;
    elseif batdate==220411
        dist_ = 1.2;
    elseif batdate==220412
        dist_ = 1.4;
    elseif batdate==220413
        dist_ = 1.2;
    elseif batdate==220414
        dist_ = 1.2;
    elseif batdate==220418
        dist_=1.4;
    else
        dist_ = 1.8;
    end
    splines_ = 14;
    [cortex_flights] = HumanBat_ClusterFlights(cortex_data,segmented_trajectories,cortex_ttl_seconds,dist_,splines_);
    %close all;  
    save(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'),'cortex_flights');
    figure(); hold on; for ii=1:length(cortex_flights.clusterIndex{2})
        numnum = cortex_flights.clusterIndex{2}(ii);
        numnum_start = cortex_flights.flight_starts_idx(numnum);
        numnum_end = cortex_flights.flight_ends_idx(numnum);
        xlim([-2.9 2.9]); ylim([-2.6 2.6]); zlim([0 2.3]);
        plot3(cortex_flights.trajectoriesContinous(1,numnum_start:numnum_end),cortex_flights.trajectoriesContinous(2,numnum_start:numnum_end),cortex_flights.trajectoriesContinous(3,numnum_start:numnum_end));
        scatter3(cortex_flights.trajectoriesContinous(1,numnum_start),cortex_flights.trajectoriesContinous(2,numnum_start),cortex_flights.trajectoriesContinous(3,numnum_start),50,'r','filled');
    end; hold off;
        %plot3(cortex_flights)
    segmented_trajectories = [];
end

%% CIHOLAS-BAT: Shift ciholas time, interpolate into 120Hz, transform to cortex 3d coordinates, segment flights
% ------------------------------------------------------------------
% Check if Cortex data has been formatted and saved
if exist(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'))
    disp("Ciholas Data has been formatted, segmented, and clustered. Loading now...");
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
    load(strcat(exp_data_path,'ciholas/clustered_ciholas_flights.mat'));
else
    disp("Segmenting Ciholas Flights")
    % Segment the ciholas flights
    [ciholas_t,ciholas_r] = HumanBat_SegmentCiholasTrajectories(ciholas_data,ciholas2cortex,cortex_data,exp_data_path,batdate);
    close all;  save(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'),'ciholas_t','ciholas_r');
    
    disp("Clustering Ciholas Flights"); clear r t r_ t
    [ciholas_flights] = HumanBat_ClusterCiholasFlights(exp_data_path,batdate);
    close all;  save(strcat(exp_data_path,'ciholas/clustered_ciholas_flights.mat'),'ciholas_flights');
end

%% CIHOLAS-HUMANS: Shift ciholas time, interpolate into 120Hz, transform to cortex 3d coordinates, segment flights
% ------------------------------------------------------------------
% Check if Cortex data has been formatted and saved
if exist(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'))
    disp("Ciholas Human Data has been formatted, segmented, and clustered. Loading now...");
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));
else
    disp("Plotting Human Trajectories")
    % posK and posM are the single-tag aligned ciholas human data, resampled to 120Hz and corrected to
    % cortex start time. Used to plot things roughly. human_t and human_r are
    % the aligned ciholas human data with all 6 human tags.
    [posK,posM,human_t,human_r,M_flights,K_flights] = HumanBat_humanPositions(ciholas_data,ciholas2cortex,cortex_data,exp_data_path,batdate);
    close all; save(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'),'human_t','human_r');
    %ciholas_data = []; cortex_data =[]; B_ephys_file =[]; cortex_file =[]; ciholas_file =[]; 
end
ciholas_data = []; cortex_data = []; B_ephys_file = []; cortex_file = []; ciholas_file = [];

% Plot Bat & KQMS positions
% figure; hold on;
% scatter(posK(:,1), posK(:,2),2,'red','filled');
% scatter(posM(:,1), posM(:,2),2,'blue','filled');
% scatter(ciholas_r(:,1), ciholas_r(:,2),2,'black','filled');
% scatter(cortex_flights.trajectoriesContinous(1,:)*1000,cortex_flights.trajectoriesContinous(2,:)*1000,2,'green','filled');
% legend('KQ', 'MS','CiholasBat','CortexBat');
% title(strcat("2d position of bats, K, and M on ",num2str(batdate)));
% hold off; 
posK = []; posM = []; ciholas_t = []; ciholas_r = []; human_t = []; human_r = []; 

%% EPHYS align to first CORTEX TTL (For now only logger 13)
% ------------------------------------------------------------------
if exist(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'))
    disp("Ephys data already aligned. Loading...");
    load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
else
    [B_ephys_data] = HumanBat_alignEphystoCortexTTL(B_ephys_data,cortex_ttl_seconds,ephys_TTL);
    save(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'),'B_ephys_data','-v7.3');
    cortex_ttl_seconds = [];
end

%% SYNCHRONY CHECK: Check if wingbeats align in the ciholas and ephys data streams
% ---------------------------------------------------------
[wingbeats_synced] = HumanBat_wingbeat_sync_check(exp_data_path,batdate,logger,ephys_TTL);
ephys_TTL = [];
if wingbeats_synced == 1
    disp("Correct Synchrony between Ciholas and Ephys!");
else
    disp("BAD SYNCHRONY! Check ciholas ephys alignment!")
end

%% SYNCHRONY CHECK #2: For Blondie-only session check if the flightpaths align between cortex and ciholas
% ----------------------------------------------------------
[cortex_ciholas_synced] = HumanBat_cortex_ciholas_sync_check();

%% For each CORTEX flight, make structure with all data streams
% ------------------------------------------------------------------
if exist(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'))
    disp("Cortex Flight Struct Exists, Loading...");
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
    load(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'));
else
    disp("Creating Cortex Flight Struct. Very Memory Consuming.");
    [cortex_flight_struct] = HumanBat_CortexFlightStructure(exp_data_path,cortex_flights,B_ephys_data,0);
    save(strcat(exp_data_path,'cortex/cortex_bat_final_flight_structure.mat'),'cortex_flight_struct','-v7.3');
    cortex_flight_struct = []; PSTH_dn = [];  PSTH_mn = [];  PSTH_up = [];  PSTH_mn = [];  PSTH_cs = []; 
end

%% For each CIHOLAS flight, make structure with all data streams
% ------------------------------------------------------------------
if exist(strcat(exp_data_path,'ciholas/ciholas_bat_final_flight_structure.mat'))
    disp("Ciholas Flight Struct Exists, Loading...");
    load(strcat(exp_data_path,'ciholas/ciholas_bat_final_flight_structure.mat'));
else
    disp("Creating Ciholas Flight Struct. Very Memory Consuming.");
    [ciholas_flight_struct] = HumanBat_CiholasFlightStructure(exp_data_path,ciholas_flights,cortex_flights,B_ephys_data,0);
    save(strcat(exp_data_path,'ciholas/ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct','-v7.3');
    ciholas_flight_struct = []; ciholas_t = [];  ciholas_r = [];  PSTH_dn = [];  PSTH_mn = [];  PSTH_up = [];  PSTH_mn = [];  PSTH_cs = []; 
end

%% For each KEVIN flight, make structure with all data streams
% Fields:
        % M position data
        % Ciholas bat 3d position
        % Cortex 3d position
        % Ephys (Ciholas bat)
        % Ephys (Cortex bat)
        % Audio
        % K position data
    % Measures:
        % Distance to K
        % Distanct to each bat
% ------------------------------------------------------------------
% [K_flight_struct] = HumanBat_KFlightStructure();
% save(strcat(exp_data_path,'ciholas/ciholas_K_final_flight_structure.mat'),'K_flight_struct');
% 

end
