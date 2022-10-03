function [SI_info] = HumanBat_SI_fx_Cortex(batdate,logger,unit,flights_to_smush)
% Function to examine the SI for each neuron (extrapolated from Nick's paper)

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

%close all;
%unit=13; 
%batdate=220404;  
%logger=15;
%flights_to_smush = 'all'; % (can also say '_1 meaning to tripod 1 or '1_' meaning from tripod 1)
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));

load(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
close all;

%% 1. Calculate 3d rate maps for flying
bin_size = 20; % cm to section room into

% a. Get room boundaries and bin into 10cmx10cmx10cm voxels
figure(); hold on; title("Room Boundaries, All Flights, No Rest"); for i=1:length(cortex_flight_struct_resort); plot3(cortex_flight_struct_resort{i}.pos(1,:),cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:)); end; hold off;
room_bounds = [-290,290; -260,260; 1,230];
room_bounds_new_origin = [room_bounds(1,:)-room_bounds(1,1);room_bounds(2,:)-room_bounds(2,1);room_bounds(3,:)];
room_2d = zeros(round(room_bounds_new_origin(1:2,2)/bin_size)');

% b. Smush all the flights together (take out the resting time) 
if strcmp(flights_to_smush,'all')
    clip = 0;
    cortex_flying_r = []; ephys_flying_r = [];
    for i=1:size(cortex_flight_struct_resort,2)
        st = cortex_flight_struct_resort{i}.fstart_trimmed+clip; ed = cortex_flight_struct_resort{i}.fend_trimmed+clip;
        cortex_flying_r = [cortex_flying_r;cortex_flights.trajectoriesContinous(:,st:ed)'];
    
        % Make binary vector of ephys for every flight
        e_bin = zeros(1,cortex_flight_struct_resort{i}.ephys_trimmed_t);
        if cortex_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
            e_bin(cortex_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
        end
        ephys_flying_r = [ephys_flying_r,e_bin];
    end
elseif contains(flights_to_smush,'_')
    if flights_to_smush(1) == '_'
        for j=1:10
            if str2double(flights_to_smush(2)) == j
                index_ = j
                disp(strcat("Concatenating all flights to location ",num2str(index_)));
            end
        end
        clip = 0;
        cortex_flying_r = []; ephys_flying_r = [];
        for i=1:size(cortex_flight_struct_resort,2)
            if cortex_flight_struct_resort{i}.tripod_landing == index_
                st = cortex_flight_struct_resort{i}.fstart_trimmed+clip; ed = cortex_flight_struct_resort{i}.fend_trimmed+clip;
                cortex_flying_r = [cortex_flying_r;cortex_flights.trajectoriesContinous(:,st:ed)];
            
                % Make binary vector of ephys for every flight
                e_bin = zeros(1,cortex_flight_struct_resort{i}.ephys_trimmed_t);
                if cortex_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
                    e_bin(cortex_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
                end
                ephys_flying_r = [ephys_flying_r,e_bin];
            end
        end
    elseif flights_to_smush(2) == '_'
        for j=1:10
            if str2double(flights_to_smush(1)) == j
                index_ = j
                disp(strcat("Concatenating all flights from location ",num2str(index_)));
            end
        end
        clip = 0;
        cortex_flying_r = []; ephys_flying_r = [];
        for i=1:size(cortex_flight_struct_resort,2)
            if cortex_flight_struct_resort{i}.tripod_takeoff == index_
                st = cortex_flight_struct_resort{i}.fstart_trimmed+clip; ed = cortex_flight_struct_resort{i}.fend_trimmed+clip;
                cortex_flying_r = [cortex_flying_r;cortex_flights.trajectoriesContinous(:,st:ed)];
            
                % Make binary vector of ephys for every flight
                e_bin = zeros(1,cortex_flight_struct_resort{i}.ephys_trimmed_t);
                if cortex_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
                    e_bin(cortex_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
                end
                ephys_flying_r = [ephys_flying_r,e_bin];
            end
        end
    end
end

% Plot concatenated flights to check
cortex_flying_r = cortex_flying_r*1000;
figure(); hold on; title(strcat(flights_to_smush," flights")); 
scatter3(cortex_flights.trajectoriesContinous(1,:),cortex_flights.trajectoriesContinous(2,:),cortex_flights.trajectoriesContinous(3,:),0.5,'k');
scatter3(cortex_flying_r(:,1),cortex_flying_r(:,2),cortex_flying_r(:,3)); hold off;

% c. Get time-spent in each voxel (occupancy map)
cortex_flying_r_rescale = [cortex_flying_r(:,1)-room_bounds(1,1)*10, cortex_flying_r(:,2)-room_bounds(2,1)*10];

% Add fake coordinates
cortex_flying_r_rescale = [cortex_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time spent in each 20x20 voxel'); hist3(cortex_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap,Xedges,Yedges,binX,binY] = histcounts2(cortex_flying_r_rescale(:,1),cortex_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
%figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
occMap_ts = occMap/120;  figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
occMap_ts_sm = imgaussfilt(occMap_ts);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
occMap_rollout = reshape(occMap_ts,[1,size(occMap,1)*size(occMap,2)]);

% d. For every XY-bin, get the # of spikes in that bin (1xn vector of
% spikes where n is the cumulative length of occupied bins)
coord_pairs = [binX,binY]; coord_pairs_unique = unique(coord_pairs,'rows');
ephys_binCount = []; 
for j=1:length(coord_pairs_unique)
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
            temp_count = temp_count + sum(ephys_flying_r((i-1)*floor(length(ephys_flying_r)/length(cortex_flying_r_rescale))+1:(i)*floor(length(ephys_flying_r)/length(cortex_flying_r_rescale))));
        end
    end
    ephys_binCount(j) = temp_count;
end
save(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit',num2str(unit),'_binCount_flights',flights_to_smush),'ephys_binCount','coord_pairs','coord_pairs_unique','occMap','binX','binY','-v7.3')

%% f. Divide the number of spikes in that XY-bin by the time spent in that bin
% ephys_binCount_norm = zeros(1,length(ephys_binCount)); occ_bin = []; occ_bin_ts=[]; 
% for i=1:length(ephys_binCount)
%     aa = find(all(coord_pairs==coord_pairs_unique(i,:),2));
%     occ_bin(i) = length(aa);
%     occ_bin_ts(i) = length(aa)/120;
%     if occ_bin_ts(i) < 0.005
%         occ_bin_ts(i)=0;
%     end
%     ephys_binCount_norm(i) = ephys_binCount(i)/occ_bin_ts(i);
% end

%% 2. Set occupancy map bins to 0 if the bat spent less than 10ms in that bin
for i=1:length(occMap_rollout)
    if occMap_rollout(i) < 0.2
        occMap_rollout(i) = 0;
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
        for m=1:length(ephys_binCount)
            if coord_pairs_unique(m,1) == ubinX(i) & coord_pairs_unique(m,2) == ubinY(j)
                %spike_map_norm_sm(i,j) = ephys_binCount_norm_sm(m);
                %spike_map_norm(i,j) = ephys_binCount_norm(m);
                spike_map_raw(ubinX(i),ubinY(j)) = ephys_binCount(m);
            end
        end
    end
end

% % 3. Visualize the 2d ratemap
% spike_map_raw = zeros(length(unique(binX)),length(unique(binY)));  
% spike_map_raw = zeros(size(occMap_ts_thresh,1),size(occMap_ts_thresh,2));
% for i=1:length(unique(binX))
%     for j=1:length(unique(binY))
%         for m=1:length(ephys_binCount)
%             if coord_pairs_unique(m,1) == i & coord_pairs_unique(m,2) == j
%                 %spike_map_norm_sm(i,j) = ephys_binCount_norm_sm(m);
%                 %spike_map_norm(i,j) = ephys_binCount_norm(m);
%                 spike_map_raw(i,j) = ephys_binCount(m);
%             end
%         end
%     end
% end

% a. Calculate the raw place map using the absolute number of spikes per bin divided by the occupancy in seconds in each bin
place_map_raw = spike_map_raw./occMap_ts_thresh;
figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap raw")); hold on; imagesc(place_map_raw'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap raw")); hold off;

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
f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap")); hold on; colorbar; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN")); hold off;
%saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap',flights_to_smush),'jpg');
SI_info.information_per_spike = information_per_spike;
SI_info.unit = unit;
SI_info.logger = logger;
SI_info.batdate = batdate;

%% 5. Check if the shuffles already exist for that cell
% If they exist, run a t test to assess significance of that unit's SI
if exist(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI.mat'))
    SI_info.shuffled = 1;
    shuffles = load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI.mat'));
    [h,p] = ttest2(information_per_spike,shuffles.information_per_spike);
    SI_info.h = h;
    SI_info.p = p;
    if h==1 & information_per_spike < mean(shuffles.information_per_spike) 
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike);
        SI_info.SI_category = 'SI < shuffle SI';
        disp("Unit has less spatial information than the mean of the shuffles.");
        figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike)+information_per_spike-max(shuffles.information_per_spike)+0.1]);
    elseif h==1 & mean(shuffles.information_per_spike) < information_per_spike
        SI_info.SI_category = 'Spatially Tuned. SI > shuffle SI';
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike);
        disp("This Unit has significant spatial information.");
        disp(strcat("SI of unit ",num2str(unit),":",num2str(information_per_spike)));
        figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike)+information_per_spike-max(shuffles.information_per_spike)+0.1]);
    else
       SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike);
       SI_info.SI_category = 'None';
       disp("Not spatially tuned unit.");
       figure();  hold on; title(strcat("Spatial information of shuffles and Unit ",num2str(unit)," (red dot)")); boxplot(shuffles.information_per_spike); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike)+information_per_spike-max(shuffles.information_per_spike)+0.1]);
    end
else
    SI_info.shuffled = 0; 
    disp("Run Shuffling on Seperate Matlab instance! (HumanBat_shuffle_SI.m)");
end

end