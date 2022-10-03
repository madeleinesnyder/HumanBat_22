%function [SI_info] = HumanBat_SI_OBPosition_fx(batdate,logger,unit)
% Function to examine the SI for each neuron given a map of the HUMAN MOVEMENT 
% (extrapolated from Nick's paper)

% Inputs:
    % unit: putative unit to examine SI of (i.e. 5)
    % batdate: date of session
    % logger: which logger data to examine
    % bin_size: the size of the bins to segment the room into (i.e. 20 (cm))
% Outputs:
    % Significance of spatial information for that unit
    % Heatmap of the unit and the shuffles

% MCS 6/22/22
% ============================================

% ALWAYS USE LOGGER 13 FOR THIS SCRIPT
logger=13;

exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));

load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'ciholas/clustered_ciholas_flights.mat'));
load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat')); 
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'))
close all;

%% 1. Calculate 3d rate maps for flying
bin_size = 20; % cm to section room into

OB = cortex_flights.trajectoriesContinous';

% a. Get room boundaries and bin into 10cmx10cmx10cm voxels
%figure(); hold on; title("Room Boundaries, All oTHER bAT Movement and Rest"); plot3(OB(:,1),OB(:,2),OB(:,3)); hold off;
room_bounds = [-290,290; -260,260; 1,230];
room_bounds_new_origin = [room_bounds(1,:)-room_bounds(1,1);room_bounds(2,:)-room_bounds(2,1);room_bounds(3,:)];
room_2d = zeros(round(room_bounds_new_origin(1:2,2)/bin_size)');

%% b. Smush all the flights together (take out the resting time) 

% Find all time that ciholas bat was chillin in a spot 
spots = [-2563,1366,1714;
         1920,-2005,1513];

bat_in_spot1_idxs = [];
for i=1:length(ciholas_r)
    if pdist2(ciholas_r(i,:),spots(1,:)) < 300
        bat_in_spot1_idxs = [bat_in_spot1_idxs,i];
    end
end
figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]); plot3(ciholas_r(bat_in_spot1_idxs,1),ciholas_r(bat_in_spot1_idxs,2),ciholas_r(bat_in_spot1_idxs,3));

bs1 = zeros(1,length(ciholas_r));
bs1(bat_in_spot1_idxs) = 1;
cs1 = zeros(1,length(ciholas_r));
crp = [];
for i=1:length(cortex_flight_struct_resort)
    crp = [crp,cortex_flight_struct_resort{i}.fstart_trimmed];
end
cs1(crp) = 1;
figure(); hold on; title("Blue is periods of ciholasbat rest, red is periods of cortex bat flight"); stem(bs1); stem(cs1);

disp(strcat("Ciholas Bat was Chillin in favorite spot for ",num2str(length(bat_in_spot1_idxs)/120/60) ," minutes over the whole session"));

% Smush all flights from cortex bat that occurred during chillin times for ciholas bat into one vector
static_cortex_flight_idxs_ = []; 
for i=1:length(cortex_flight_struct_resort)
    flight_timevec = [cortex_flight_struct_resort{i}.fstart_trimmed:cortex_flight_struct_resort{i}.fend_trimmed];
    static_flight = 0;
    for j=1:length(flight_timevec)
        if ismember(flight_timevec(j),bat_in_spot1_idxs)
            static_flight = static_flight+1;
        else 
            static_flight = static_flight;
        end
    end
    if static_flight == 0
        disp("this flight didn't occur during a period of stillness for ciholas bat")
    elseif static_flight == length(flight_timevec)
        static_cortex_flight_idxs_ = [static_cortex_flight_idxs_,i];
    elseif static_flight < length(flight_timevec)-60
        disp("partial stillnesss during this cortex flight. plot")
    end
end

% Combine all flights of cortexbat that occurred during ciholasbat
% stillness
OB_flying_r = []; OB_flying_all_r = []; lgc=0; static_cortex_flight_idxs = []; 
for i=1:length(static_cortex_flight_idxs_)
    %if cortex_flight_struct_resort{static_cortex_flight_idxs_(i)}.fclus ~=1
        OB_flying_r = [OB_flying_r, cortex_flight_struct_resort{static_cortex_flight_idxs_(i)}.pos];
        lgc = lgc+1;
        static_cortex_flight_idxs = [static_cortex_flight_idxs,static_cortex_flight_idxs_(i)];
    %else
    %    OB_flying_all_r = [OB_flying_all_r, cortex_flight_struct_resort{static_cortex_flight_idxs_(i)}.pos];
    %end
end
OB_flying_r = OB_flying_r';

disp(strcat("During the CiholasBat's X minutes of rest, there were ", num2str(size(OB_flying_r,1)/120)," seconds of CLUSTERABLE flying by the CortexBat, aka ",num2str(lgc)," flights."));

%% Extract the periods of ephys for CiholasBat from the start and stop times
% of CortexBat
ephys_flying_r = []; ciholas_ephys_on_cortex = {};
for i=1:length(static_cortex_flight_idxs)
    flight_num = static_cortex_flight_idxs(i);
    start_time = cortex_flight_struct_resort{flight_num}.fstart_trimmed;
    end_time = cortex_flight_struct_resort{flight_num}.fend_trimmed;
    ephys_f_start = start_time/120*1e6;
    ephys_f_end = end_time/120*1e6;
    
    % Make a vector tt that is the length of the seconds 1e6 in that flight
    tt{flight_num} = [0:(ephys_f_end-ephys_f_start)];   
    temp_rr = zeros(length(tt{flight_num}),1);
    temp_rr_length_check(i) = length(temp_rr);

    % Fill in the rest of the timestamps with ephys_f_start
    % subtracted for that ephys stretch
    [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
    if isempty(ephys_idxs_in_range)
        disp("This unit has no activity during this CortexBat flight");
        ephys_flying_r = [ephys_flying_r;temp_rr]; %ciholas_flight_struct{flight_num}.ephys{i} = 0;
        ciholas_ephys_on_cortex{flight_num} = temp_rr;
    else
        disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during flight ",num2str(flight_num)));
        for j=1:length(ephys_idxs_in_range)
            temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
            temp_rr(round(temp_ts-ephys_f_start)) = 1;
        end
        ciholas_ephys_on_cortex{flight_num} = temp_rr;
        ephys_flying_r = [ephys_flying_r;temp_rr];
    end 
end

ciholas_flying_r = ciholas_r(bat_in_spot1_idxs,:);

% Check in on where the bat went during this time
figure(); hold on; title("Red is ciholas bat location, blue is cortexbat location");
scatter3(OB_flying_r(:,1),OB_flying_r(:,2),OB_flying_r(:,3),3,'b','filled');
plot3(ciholas_flying_r(:,1),ciholas_flying_r(:,2),ciholas_flying_r(:,3),'k');

% c. Get time-spent in each voxel (occupancy map)
OB_flying_r_rescale = [OB_flying_r(:,1)-room_bounds(1,1)*10, OB_flying_r(:,2)-room_bounds(2,1)*10];

% Add fake coordinates
OB_flying_r_rescale = [OB_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time M spent in each 20x20 voxel'); hist3(OB_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap,Xedges,Yedges,binX,binY] = histcounts2(OB_flying_r_rescale(:,1),OB_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
%figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
occMap_ts = occMap/120;  figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts)); 
occMap_ts_sm = imgaussfilt(occMap_ts);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
occMap_rollout = reshape(occMap_ts,[1,size(occMap,1)*size(occMap,2)]);

%% d. For every XY-bin, get the # of spikes in that bin (1xn vector of
% spikes where n is the cumulative length of occupied bins)
coord_pairs = [binX,binY]; coord_pairs_unique = unique(coord_pairs,'rows');
ephys_binCount = []; 
for j=1:length(coord_pairs_unique)
%     if coord_pairs_unique(j,:) == [17 21]
%         disp("The weird hotspot")
%     end
    if j == round(length(coord_pairs_unique)/4)
        disp("25% Finished");
    elseif j == round(length(coord_pairs_unique)/2)
        disp("50% Finished");
    elseif j == round(length(coord_pairs_unique)/4*3)
        disp("75% Finished");
    end
    temp_count = 0;
    for i=1:length(coord_pairs)-1
        if coord_pairs(i,:) == coord_pairs_unique(j,:)
            temp_count = temp_count + sum(ephys_flying_r((i-1)*floor(length(ephys_flying_r)/length(OB_flying_r_rescale))+1:(i)*floor(length(ephys_flying_r)/length(OB_flying_r_rescale))));
        end
    end
    ephys_binCount(j) = temp_count;
end
save(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit',num2str(unit),'_binCount_OB_traverse'),'ephys_binCount','coord_pairs','coord_pairs_unique','occMap','binX','binY','-v7.3')

%% 2. Set occupancy map bins to 0 if the bat spent less than 10ms in that bin
for i=1:size(occMap_ts,1)
    for j=1:size(occMap_ts,2)
        if occMap_ts(i,j) < 0.2
            occMap_ts(i,j) = 0;
        end
    end
end
occMap_ts_thresh = reshape(occMap_rollout,[size(occMap,1),size(occMap,2)]); 
occMap_ts_thresh = occMap_ts; 

% 3. Visualize the 2d ratemap
spike_map_raw = zeros(length(unique(binX)),length(unique(binY)));  
spike_map_raw = zeros(size(occMap_ts_thresh,1),size(occMap_ts_thresh,2));
ubinX = unique(binX);
ubinY = unique(binY);
for i=1:length(unique(binX))
    for j=1:length(unique(binY))
%         if i==21 & j==16
%             disp("Hot spot")
%         end
        for m=1:length(ephys_binCount)
            if coord_pairs_unique(m,1) == ubinX(i) & coord_pairs_unique(m,2) == ubinY(j)
                %spike_map_norm_sm(i,j) = ephys_binCount_norm_sm(m);
                %spike_map_norm(i,j) = ephys_binCount_norm(m);
                spike_map_raw(ubinX(i),ubinY(j)) = ephys_binCount(m);
            end
        end
    end
end

% a. Calculate the raw place map using the absolute number of spikes per bin divided by the occupancy in seconds in each bin

place_map_raw = spike_map_raw./occMap_ts_thresh;
figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap raw of OB trajectories while this brain is stationary")); hold on; imagesc(place_map_raw'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap of stationary brain while other bat flying; raw")); hold off;
figure(); hold on; imagesc(place_map_raw'); axis off; hold off;

% b. Create gaussian kernel to smooth the raw spike map
sigma = 1.5; % Hafting I guess?
hsize = 5*round(sigma)+1;
gaussian_kernel = fspecial('gaussian',hsize,sigma);
%figure(); surf(gaussian_kernel); % If you want to see how your kernel looks like 
spike_map_smoothed = imfilter(spike_map_raw,gaussian_kernel);
occMap_ts_thresh_smoothed = imfilter(occMap_ts_thresh,gaussian_kernel);
place_map_smoothed = spike_map_smoothed./ occMap_ts_thresh_smoothed ;
figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap smoothed with gaussian kernel")); hold on; imagesc((place_map_smoothed')); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smoothed with gaussian kernel"))

% c. Account for nearest neighbor activations in a new place map
% From Michael's hippo code
idx_timespent_density = zeros(size(occMap,1),size(occMap,2)) + NaN ; % Initialize
for ii_y_bin = 2 : size(occMap,1)-1 %length(unique(binX)) - 1 , % Loop over x-bins
    for ii_x_bin = 2 : size(occMap,2)-1 %length(unique(binY)) - 1 , % Loop over y-bins
        matrix_3x3_of_neighbors = ...
            occMap_ts( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
        sum_including_the_central_bin = sum(sum( matrix_3x3_of_neighbors)); % Count the matrix_3x3_of_neighbors + the central bin itself
        if ( sum_including_the_central_bin  > 0 ), % If the animal visited any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_y_bin,ii_x_bin) = 1; % Put 1 in the central bin
        else  % If the animal did NOT visit any of this bin's 8 neighbors (3x3 region)
            idx_timespent_density(ii_y_bin,ii_x_bin) = 0; % Put 0 in the central bin (later we will divide by this 0 and will get NaN for the firing-rate map)
        end
    end
end

% Smooth the place map with NaNs
place_map_smoothed_with_NaN = (place_map_smoothed.* idx_timespent_density)./idx_timespent_density;
occMap_ts_thresh_smoothed_with_NaN = (occMap_ts_thresh_smoothed.* idx_timespent_density)./idx_timespent_density;

idx_notNaN_PlaceField = find( ~isnan( place_map_smoothed_with_NaN ) ); % Find the indexes of non-NaN bins
idx_isNaN_PlaceField = find( isnan( place_map_smoothed_with_NaN  ) ); % Find the indexes of NaN bins

idx_notNaN_PlaceField_un_smoothed_rate_map = find( ~isnan( place_map_raw ) ); % Find the indexes of non-NaN bins
figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap")); hold on; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN")); hold off;

% Janky fix:
% for i=1:size(place_map_smoothed_with_NaN,1)
%     for j=1:size(place_map_smoothed_with_NaN,2)
%         if place_map_smoothed_with_NaN(i,j) > nanmean(nanmean(place_map_smoothed_with_NaN))*2
%             place_map_smoothed_with_NaN(i,j) = nanmean(nanmean(place_map_smoothed_with_NaN));
%         end
%     end
% end
% figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap with outliers removed")); hold on; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN")); hold off;


%% 4. Calculate INFORMATION PER SPIKE  (computed for the SMOOTHED field): -----
% Information_per_spike = sum( p_i * ( r_i / r ) * log2( r_i / r ) )
%    Where:
%       r_i = firing rate in bin i ;
%       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
%       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
% See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

r_i = place_map_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED Place Field
p_i = occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
    sum( occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
% % %             r = mean( r_i ) ;
r = sum( r_i .* p_i );
information_per_spike = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)

% Save figure of the heatmap 
f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap for OB position")); hold on; colorbar; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN for OB position while this brain is stationary")); hold off;

%saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap',flights_to_smush),'jpg');

%% Plot the red dots for this particular case 

figure(); hold on; title(strcat("Flights of other bat with ephys plotted on top")); xlim([-2900 2900]); ylim([-2500 2500]); zlim([0 2300]);

for pp=1:length(static_cortex_flight_idxs)
    i = static_cortex_flight_idxs(pp);

    scatter3(cortex_flight_struct_resort{i}.pos(1,:),cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:),2,'b','filled');

    % Transform ephys timestamps to ciholas time
    % Make artificial vector the time duration of ciholas
    flight_dur = (cortex_flight_struct_resort{i}.fend_trimmed - cortex_flight_struct_resort{i}.fstart_trimmed);
    flight_vec = [1:flight_dur];
    if sum(ciholas_ephys_on_cortex{i}) == 0
        disp(strcat("No unit firing on flight ",num2str(i)));
        %continue;
    else
        ciholas_ephys_timestamps_in_ciholas_time = find(ciholas_ephys_on_cortex{i}==1)/1e6*120;
        for j=1:length(ciholas_ephys_timestamps_in_ciholas_time)
            nearest_ephys_point = dsearchn(flight_vec',ciholas_ephys_timestamps_in_ciholas_time(j));
            scatter3(cortex_flight_struct_resort{i}.pos(1,nearest_ephys_point),cortex_flight_struct_resort{i}.pos(2,nearest_ephys_point),cortex_flight_struct_resort{i}.pos(3,nearest_ephys_point),20,'r','filled');          
        end
    end
end

% Checks
figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
for i=1:length(static_cortex_flight_idxs)
    stt = cortex_flight_struct_resort{static_cortex_flight_idxs(i)}.fstart_trimmed;
    ett = cortex_flight_struct_resort{static_cortex_flight_idxs(i)}.fend_trimmed;
    plot3(ciholas_r(stt:ett,1),ciholas_r(stt:ett,2),ciholas_r(stt:ett,3),'b');
    plot3(cortex_flight_struct_resort{static_cortex_flight_idxs(i)}.pos(1,:),cortex_flight_struct_resort{static_cortex_flight_idxs(i)}.pos(2,:),cortex_flight_struct_resort{static_cortex_flight_idxs(i)}.pos(3,:),'r');
end

counter = 0; for kk = 1:length(ciholas_ephys_on_cortex); counter = counter + sum(ciholas_ephys_on_cortex{kk}); end;



%% 5. Check if the shuffles already exist for that cell
% If they exist, run a t test to assess significance of that unit's SI
if exist(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit_',num2str(unit),'_shuffles_OB/shuffle_SI_OB.mat'))
    shuffles = load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit_',num2str(unit),'_shuffles_OB/shuffle_SI_OB.mat'));
    [h,p] = ttest2(information_per_spike,shuffles.information_per_spike_shuff);
    if h==1 & information_per_spike < mean(shuffles.information_per_spike_shuff) 
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
        SI_info.SI_category = 'SI < shuffle SI';
        disp("Unit has less spatial information than the mean of the shuffles.");
        figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
    elseif h==1 & mean(shuffles.information_per_spike_shuff) < information_per_spike
        SI_info.SI_category = 'Spatially Tuned. SI > shuffle SI';
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
        disp("This Unit has significant spatial information.");
        disp(strcat("SI of unit ",num2str(unit),":",num2str(information_per_spike)));
        figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
    else
       SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
       SI_info.SI_category = 'None';
       disp("Not spatially tuned unit.");
       figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
    end
else
    disp("Running Shuffling!");

    num_shuffles=500;

    % Clock rotation shuffle
    % Generate all the shuffle moves:
    shuffle_clock_times = round(rand(num_shuffles,1),5);
    shuffle_clock_times = unique(shuffle_clock_times);
    for lb = 1:length(shuffle_clock_times)
        disp(strcat("Shuffle #",num2str(lb)));
        % Pick clock direction
        cs_amount = length(ephys_flying_r)*shuffle_clock_times(lb);
        round(cs_amount);
        ephys_flying_r_shift = circshift(ephys_flying_r,round(cs_amount));
        ephys_binCount_shuff = []; 
        for j=1:length(coord_pairs_unique)
            %disp(j);
            if j == round(length(coord_pairs_unique)/4)
                disp("25% Finished");
            elseif j == round(length(coord_pairs_unique)/2)
                disp("50% Finished");
            elseif j == round(length(coord_pairs_unique)/4*3)
                disp("75% Finished");
            end
            temp_count = 0;
            for i=1:length(coord_pairs)-3
                if coord_pairs(i,:) == coord_pairs_unique(j,:)
                    temp_count = temp_count + sum(ephys_flying_r_shift((i-1)*round(length(ephys_flying_r_shift)/length(OB_flying_r_rescale))+1:(i)*round(length(ephys_flying_r_shift)/length(OB_flying_r_rescale))));
                end
            end
            ephys_binCount_shuff(j) = temp_count;
        end
    
        %% Set occupancy map bins to 0 if the bat spent less than 10ms in that bin
        for i=1:length(occMap_rollout)
            if occMap_rollout(i) < 0.2
                occMap_rollout(i) = 0;
            end
        end

        %occMap_ts_thresh = reshape(occMap_rollout,[size(occMap,1),size(occMap,2)]);
        occMap_ts_thresh = occMap_ts;

        % Visualize the 2d ratemap
        spike_map_raw = zeros(length(unique(binX)),length(unique(binY)));  
        spike_map_raw = zeros(size(occMap_ts_thresh,1),size(occMap_ts_thresh,2));
        ubinX = unique(binX);
        ubinY = unique(binY);
        for i=1:length(unique(binX))
            for j=1:length(unique(binY))
                for m=1:length(ephys_binCount_shuff)
                    if coord_pairs_unique(m,1) == ubinX(i) & coord_pairs_unique(m,2) == ubinY(j)
                        spike_map_raw(ubinX(i),ubinY(j)) = ephys_binCount_shuff(m);
                    end
                end
            end
        end
        
        % Calculate the raw place map using the absolute number of spikes per bin divided by
        % the occupancy in seconds in each bin
        clear place_map_raw;
        place_map_raw = spike_map_raw./occMap_ts_thresh;
        
        % Create gaussian kernel: 
        sigma = 1.5; % Hafting I guess?
        hsize = 5*round(sigma)+1;
        gaussian_kernel = fspecial('gaussian',hsize,sigma);
        
        % Smoothing = convolve with gaussian kernel: 
        spike_map_smoothed = imfilter(spike_map_raw,gaussian_kernel);
        occMap_ts_thresh_smoothed = imfilter(occMap_ts_thresh,gaussian_kernel);
        place_map_smoothed = spike_map_smoothed./ occMap_ts_thresh_smoothed ;
        
        % Account for nearest neighbor activations in a new place map
        idx_timespent_density = zeros(size(occMap,1),size(occMap,2)) + NaN ; % Initialize
        for ii_y_bin = 2 : length(unique(binX)) - 1 , % Loop over x-bins
            for ii_x_bin = 2 : length(unique(binY)) - 1 , % Loop over y-bins
                matrix_3x3_of_neighbors = ...
                    occMap_ts( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
                sum_including_the_central_bin = sum(sum( matrix_3x3_of_neighbors)); % Count the matrix_3x3_of_neighbors + the central bin itself
                if ( sum_including_the_central_bin  > 0 ), % If the animal visited any of this bin's 8 neighbors (3x3 region)
                    idx_timespent_density(ii_y_bin,ii_x_bin) = 1; % Put 1 in the central bin
                else  % If the animal did NOT visit any of this bin's 8 neighbors (3x3 region)
                    idx_timespent_density(ii_y_bin,ii_x_bin) = 0; % Put 0 in the central bin (later we will divide by this 0 and will get NaN for the firing-rate map)
                end
            end
        end
        
        clear place_map_smoothed_with_NaN occMap_ts_thresh_smoothed_with_NaN
        place_map_smoothed_with_NaN = (place_map_smoothed.* idx_timespent_density)./idx_timespent_density;
        occMap_ts_thresh_smoothed_with_NaN = (occMap_ts_thresh_smoothed.* idx_timespent_density)./idx_timespent_density;
        
        idx_notNaN_PlaceField = find( ~isnan( place_map_smoothed_with_NaN ) ); % Find the indexes of non-NaN bins
        idx_isNaN_PlaceField = find( isnan( place_map_smoothed_with_NaN  ) ); % Find the indexes of NaN bins
        
        idx_notNaN_PlaceField_un_smoothed_rate_map = find( ~isnan( place_map_raw ) ); % Find the indexes of non-NaN bins

        % Compute SI for this shuffle
        clear r_i p_i r
        r_i = place_map_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED Place Field
        p_i = occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
            sum( occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
        r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
        % % %             r = mean( r_i ) ;
        r = sum( r_i .* p_i );
        information_per_spike_shuff(lb) = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)
    
        f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap shuffle ",num2str(lb)," I: ",num2str(information_per_spike_shuff(lb))),'visible','off'); hold on; imagesc(place_map_smoothed_with_NaN'); title(strcat("Unit ",num2str(unit)," 2D ratemap shuffle OB ",num2str(lb)," I: ",num2str(information_per_spike_shuff(lb)))); hold off;
        if ~exist(strcat(exp_data_path,"ephys/logger13/extracted_data/unit_",num2str(unit),"_shuffles_OB"))
            mkdir(strcat(exp_data_path,"ephys/logger13/extracted_data/unit_",num2str(unit),"_shuffles_OB"));
            saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles_OB/',num2str(batdate),'_14650_PlaceMap_shuffle_OB_',num2str(lb)),'jpg');
        else
            saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles_OB/',num2str(batdate),'_14650_PlaceMap_shuffle_OB_',num2str(lb)),'jpg');
        end
    end
    save(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles_OB/shuffle_SI_OB.mat'),'information_per_spike_shuff','-v7.3')
    
    disp("Now doing sig test with new shuffled data!");
    shuffles = load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit_',num2str(unit),'_shuffles_OB/shuffle_SI_OB.mat'));
    [h,p] = ttest2(information_per_spike,shuffles.information_per_spike_shuff);
    if h==1 & information_per_spike < mean(shuffles.information_per_spike_shuff) 
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
        SI_info.SI_category = 'SI < shuffle SI';
        disp("Unit has less spatial information than the mean of the shuffles.");
        figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
    elseif h==1 & mean(shuffles.information_per_spike_shuff) < information_per_spike
        SI_info.SI_category = 'Spatially Tuned. SI > shuffle SI';
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
        disp("This Unit has significant spatial information.");
        disp(strcat("SI of unit ",num2str(unit),":",num2str(information_per_spike)));
        figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
    else
       SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
       SI_info.SI_category = 'None';
       disp("Not spatially tuned unit.");
       figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
    end
end


