function [] = HumanBat_SI_HumanPosition_SharedPixels()

exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));

load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'ciholas/clustered_ciholas_flights.mat'));
load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat')); 
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'))
close all;

%% 1. Calculate 3d rate maps for flying
bin_size = 20; % cm to section room into

% a. Get room boundaries and bin into 10cmx10cmx10cm voxels
figure(); hold on; title("Room Boundaries, All Human Movement and Rest"); plot3(human_r(:,1,3),human_r(:,2,3),human_r(:,3,3),'r'); plot3(human_r(:,1,4),human_r(:,2,4),human_r(:,3,4),'b'); hold off;
room_bounds = [-290,290; -260,260; 1,230];
room_bounds_new_origin = [room_bounds(1,:)-room_bounds(1,1);room_bounds(2,:)-room_bounds(2,1);room_bounds(3,:)];
room_2d = zeros(round(room_bounds_new_origin(1:2,2)/bin_size)');

%% b. Smush all the flights together (take out the resting time) 

% Smushing only human approach bat time
human_flying_r_ = human_r; 
if B_ephys_data.TT_unit(unit).AlignedTimestamps(1) < 0
    ephys_flying_maximum = B_ephys_data.TT_unit(unit).AlignedTimestamps(end)-B_ephys_data.TT_unit(unit).AlignedTimestamps(1);
else
    ephys_flying_maximum = round(B_ephys_data.TT_unit(unit).AlignedTimestamps(end));
end

ephys_flying_r_ = zeros(1,ephys_flying_maximum);
for i=1:length(B_ephys_data.TT_unit(unit).AlignedTimestamps)
    if B_ephys_data.TT_unit(unit).AlignedTimestamps(i) > 0
        ephys_flying_r_(round(B_ephys_data.TT_unit(unit).AlignedTimestamps(i))) = 1;
    end
end

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

%% Make an occupancy map for each human
madeleine_flying_r = human_flying_r_(bat_in_spot1_idxs,:,4);
kevin_flying_r = human_flying_r_(bat_in_spot1_idxs,:,3);

figure(); hold on; plot3(madeleine_flying_r(:,1),madeleine_flying_r(:,2),madeleine_flying_r(:,3),'b'); plot3(kevin_flying_r(:,1),kevin_flying_r(:,2),kevin_flying_r(:,3),'r');

madeleine_flying_r_rescale = [madeleine_flying_r(:,1)-room_bounds(1,1)*10, madeleine_flying_r(:,2)-room_bounds(2,1)*10];
% Add fake coordinates
madeleine_flying_r_rescale = [madeleine_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time M spent in each 20x20 voxel'); hist3(madeleine_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap_M,Xedges,Yedges,binX_M,binY_M] = histcounts2(madeleine_flying_r_rescale(:,1),madeleine_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
%figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
occMap_ts_M = occMap_M/120;  figure('Name','M occupancy map'); imagesc(rot90(occMap_ts_M));
occMap_ts_sm_M = imgaussfilt(occMap_ts_M);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
occMap_rollout_M = reshape(occMap_ts_M,1,size(occMap_M,1)*size(occMap_M,2));

kevin_flying_r_rescale = [kevin_flying_r(:,1)-room_bounds(1,1)*10, kevin_flying_r(:,2)-room_bounds(2,1)*10];
% Add fake coordinates
kevin_flying_r_rescale = [kevin_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time K spent in each 20x20 voxel'); hist3(kevin_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap_K,Xedges,Yedges,binX_K,binY_K] = histcounts2(kevin_flying_r_rescale(:,1),kevin_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
%figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
occMap_ts_K = occMap_K/120;  figure('Name','K occupancy map'); imagesc(rot90(occMap_ts_K));
occMap_ts_sm_K = imgaussfilt(occMap_ts_K);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
occMap_rollout_K = reshape(occMap_ts_K,1,size(occMap_K,1)*size(occMap_K,2));

% Identify the pixels  that the maps share
occMap_M_Mask = occMap_M; occMap_K_Mask = occMap_K;
occMap_M_Mask(occMap_M > 0) = 1;   occMap_K_Mask(occMap_K > 0) = 1; 
occMap_M_Mask_ro = reshape(occMap_M_Mask,1,size(occMap_M_Mask,1)*size(occMap_M_Mask,2));
occMap_K_Mask_ro = reshape(occMap_K_Mask,1,size(occMap_K_Mask,1)*size(occMap_K_Mask,2));
occMap_shared_mask = zeros(1,size(occMap_K_Mask_ro,2));

for i=1:length(occMap_K_Mask_ro)
    if occMap_M_Mask_ro(i) == 1 & occMap_K_Mask_ro(i) == 1
        occMap_shared_mask(i) = 2;
    else
        occMap_shared_mask(i) = 0;
    end
end

for i=1:length(occMap_shared_mask)
    if occMap_shared_mask(i) == 2
        occMap_M_shared(i) = occMap_rollout_M(i);
    else
        occMap_M_shared(i) = 0;
    end
end

for i=1:length(occMap_shared_mask)
    if occMap_shared_mask(i) == 2
        occMap_K_shared(i) = occMap_rollout_K(i);
    else
        occMap_K_shared(i) = 0;
    end
end

% Occupancy maps for each human, using only the shared pixels (just occMap
% but pruned)
occMap_M_sh = reshape(occMap_M_shared,size(occMap_M_Mask,1),size(occMap_M_Mask,2));
occMap_K_sh = reshape(occMap_K_shared,size(occMap_K_Mask,1),size(occMap_K_Mask,2));


%% For each human, find the set of cortex timepoints where the occupancy was shared

sm = reshape(occMap_shared_mask,size(occMap_M_Mask,1),size(occMap_M_Mask,2)); sm(sm==2)=1;
coord_pairs_unique = [];
for i=1:size(sm,1)
    for j=1:size(sm,2)
        if sm(i,j) == 1
            coord_pairs_unique = [coord_pairs_unique;[i j]];
        end
    end
end
% Find timepoints in M that corrospond to being in shared pixels
tp_include_M = [];
coord_pairs = [binX_M,binY_M]; 
for j=1:length(coord_pairs_unique)
    for i=1:length(coord_pairs)-1
        if coord_pairs(i,:) == coord_pairs_unique(j,:)
            tp_include_M = [tp_include_M,i];
        end
    end
end
% Find timepoints in K that corrospond to being in shared pixels
clear coord_pairs 
tp_include_K = [];
coord_pairs = [binX_K,binY_K]; 
for j=1:length(coord_pairs_unique)
    for i=1:length(coord_pairs)-1
        if coord_pairs(i,:) == coord_pairs_unique(j,:)
            tp_include_K = [tp_include_K,i];
        end
    end
end
bhh = [];
for i=1:length(bat_in_spot1_idxs)
    if ismember(bat_in_spot1_idxs(i),tp_include_K) & ismember(bat_in_spot1_idxs(i),tp_include_M)
        bhh = [bhh,bat_in_spot1_idxs(i)];
    end
end

% Plot the raw human data with these timepoints to see what's what
figure(); hold on;
scatter3(human_r(bhh,1,4),human_r(bhh,2,4),human_r(bhh,3,4),'b');
scatter3(human_r(bhh,1,3),human_r(bhh,2,3),human_r(bhh,3,3),'r');
scatter3(ciholas_r(bhh,1),ciholas_r(bhh,2),ciholas_r(bhh,3),'k');

%% Extract ephys data using segmenting method
in_stretch = 0; stretch = []; stretches = {}; individuals = {};
for i=1:length(bhh)-1
    if bhh(i+1)-bhh(i) == 1
        in_stretch = 1;
    else
        in_stretch = 0;
    end
    if in_stretch == 1 & ~isempty(stretch)
        stretch = [stretch,bhh(i)];
    elseif in_stretch == 1 & isempty(stretch)
        stretch = [stretch, bhh(i)];
    elseif in_stretch == 0 & ~isempty(stretch)
        stretches{end+1} = [stretch];
        stretch = [];
    elseif in_stretch == 0 & isempty(stretch)
        individuals{end+1} = bhh(i);
        stretch = [];
    end
end

% Check if stretches add up
running_tally = 0;
for i=1:length(stretches)
    running_tally = running_tally+length(stretches{i});
end

% For each stretch, make the ephys!
ephys_flying_r = [];
for i=1:length(stretches)
    start_ = stretches{i}(1);
    end_ = stretches{i}(end);
    ephys_f_start = start_/120*1e6;
    ephys_f_end = end_/120*1e6;
    
    % Make a vector tt that is the length of the seconds 1e6 in that flight
    tt{i} = [0:(ephys_f_end-ephys_f_start)];   
    temp_rr = zeros(length(tt{i}),1);
    temp_rr_length_check(i) = length(temp_rr);

    % Fill in the rest of the timestamps with ephys_f_start
    % subtracted for that ephys stretch
    [B,ephys_idxs_in_range] = find(ephys_f_start < B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps < ephys_f_end);
    if isempty(ephys_idxs_in_range)
        disp("This unit has no activity during this stretch");
        ephys_flying_r = [ephys_flying_r;temp_rr]; %ciholas_flight_struct{flight_num}.ephys{i} = 0;
        ciholas_ephys_on_human{i} = temp_rr;
    else
        disp(strcat("This unit has ",num2str(length(ephys_idxs_in_range))," spikes during stretch ",num2str(i)));
        for j=1:length(ephys_idxs_in_range)
            temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
            temp_rr(round(temp_ts-ephys_f_start)) = 1;
        end
        ciholas_ephys_on_cortex{i} = temp_rr;
        ephys_flying_r = [ephys_flying_r;temp_rr];
    end 
end

bsp1 = [];
for i=1:length(stretches)
    bsp1 = [bsp1,stretches{i}];
end

madeleine_flying_r_bsp = madeleine_flying_r(bsp1,:);
kevin_flying_r_bsp = kevin_flying_r(bsp1,:);

figure(); hold on;
    scatter3(madeleine_flying_r_bsp(:,1),madeleine_flying_r_bsp(:,2),madeleine_flying_r_bsp(:,3),3,'b');
    scatter3(kevin_flying_r_bsp(:,1),kevin_flying_r_bsp(:,2),kevin_flying_r_bsp(:,3),3,'r');
    scatter3(ciholas_r(bsp1,1),ciholas_r(bsp1,2),ciholas_r(bsp1,3),'k');


% Now ephys_fly_r and madeliene/kevin_flying_r should be vectors of the stretches of
% continuous stillness captured by stretches{} (second half)

%% For Madeleine:
% c. Get time-spent in each voxel (occupancy map)
madeleine_flying_r_rescale = [madeleine_flying_r_bsp(:,1)-room_bounds(1,1)*10, madeleine_flying_r_bsp(:,2)-room_bounds(2,1)*10];

% Add fake coordinates
madeleine_flying_r_rescale = [madeleine_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time M spent in each 20x20 voxel'); hist3(madeleine_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap,Xedges,Yedges,binX,binY] = histcounts2(madeleine_flying_r_rescale(:,1),madeleine_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
%figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
occMap_ts = occMap/120;  figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
occMap_ts_sm = imgaussfilt(occMap_ts);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
occMap_rollout = reshape(occMap_ts,[1,size(occMap,1)*size(occMap,2)]);

%% d. For every XY-bin, get the # of spikes in that bin (1xn vector of
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
            temp_count = temp_count + sum(ephys_flying_r((i-1)*floor(length(ephys_flying_r)/length(madeleine_flying_r_rescale))+1:(i)*floor(length(ephys_flying_r)/length(madeleine_flying_r_rescale))));
        end
    end
    ephys_binCount(j) = temp_count;
end
%save(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit',num2str(unit),'_binCount_',human_name,'_traverse'),'ephys_binCount','coord_pairs','coord_pairs_unique','occMap','binX','binY','-v7.3')

%% 2. Set occupancy map bins to 0 if the bat spent less than 10ms in that bin
for i=1:length(occMap_rollout)
    if occMap_rollout(i) < 1
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
                spike_map_raw(ubinX(i),ubinY(j)) = ephys_binCount(m);
            end
        end
    end
end

% a. Calculate the raw place map using the absolute number of spikes per bin divided by the occupancy in seconds in each bin
place_map_raw = spike_map_raw./occMap_ts_thresh;
figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap raw Madeleine only shared pix")); hold on; imagesc(place_map_raw'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap raw")); hold off;

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
f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap for Madeleine position only shared pixels")); hold on; colorbar; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN for M position shared pix")); hold off;

%% Red Dots Madeleine
figure(); hold on; title("Madeleine flights with bat ephys plotted on top");
scatter3(madeleine_flying_r_bsp(:,1),madeleine_flying_r_bsp(:,2),madeleine_flying_r_bsp(:,3),3,[0.58 0.81 0.98]);

% Transform ephys timestamps to ciholas time
% Make artificial vector the time duration of ciholas
spike_idxs = find(ephys_flying_r==1)/1e6*120; % These are the indexes at which spikes occur

% For every spike, plot it at that index
for i=1:length(spike_idxs)
    scatter3(madeleine_flying_r_bsp(round(spike_idxs(i)),1),madeleine_flying_r_bsp(round(spike_idxs(i)),2),madeleine_flying_r_bsp(round(spike_idxs(i)),3),30,'r','filled');
end

%% For Kevin:
% c. Get time-spent in each voxel (occupancy map)
kevin_flying_r_rescale = [kevin_flying_r_bsp(:,1)-room_bounds(1,1)*10, kevin_flying_r_bsp(:,2)-room_bounds(2,1)*10];

% Add fake coordinates
kevin_flying_r_rescale = [kevin_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time K spent in each 20x20 voxel'); hist3(kevin_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap,Xedges,Yedges,binX,binY] = histcounts2(kevin_flying_r_rescale(:,1),kevin_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
%figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
occMap_ts = occMap/120;  figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
occMap_ts_sm = imgaussfilt(occMap_ts);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
occMap_rollout = reshape(occMap_ts,[1,size(occMap,1)*size(occMap,2)]);

%% d. For every XY-bin, get the # of spikes in that bin (1xn vector of
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
            temp_count = temp_count + sum(ephys_flying_r((i-1)*floor(length(ephys_flying_r)/length(kevin_flying_r_rescale))+1:(i)*floor(length(ephys_flying_r)/length(kevin_flying_r_rescale))));
        end
    end
    ephys_binCount(j) = temp_count;
end
%save(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/unit',num2str(unit),'_binCount_',human_name,'_traverse'),'ephys_binCount','coord_pairs','coord_pairs_unique','occMap','binX','binY','-v7.3')

%% 2. Set occupancy map bins to 0 if the bat spent less than 10ms in that bin
for i=1:length(occMap_rollout)
    if occMap_rollout(i) < 1
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
                spike_map_raw(ubinX(i),ubinY(j)) = ephys_binCount(m);
            end
        end
    end
end

% a. Calculate the raw place map using the absolute number of spikes per bin divided by the occupancy in seconds in each bin
place_map_raw = spike_map_raw./occMap_ts_thresh;
figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap raw Kevin only shared pix")); hold on; imagesc(place_map_raw'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap raw")); hold off;

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
f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap for Kevin position only shared pixels")); hold on; colorbar; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN for K position shared pix")); hold off;

%% Red Dots Kevin
figure(); hold on; title("Kevin flights with bat ephys plotted on top");
scatter3(kevin_flying_r_bsp(:,1),kevin_flying_r_bsp(:,2),kevin_flying_r_bsp(:,3),3,[1 0.6 0.9]);

% Transform ephys timestamps to ciholas time
% Make artificial vector the time duration of ciholas
spike_idxs = find(ephys_flying_r==1)/1e6*120; % These are the indexes at which spikes occur

% For every spike, plot it at that index
for i=1:length(spike_idxs)
    scatter3(kevin_flying_r_bsp(round(spike_idxs(i)),1),kevin_flying_r_bsp(round(spike_idxs(i)),2),kevin_flying_r_bsp(round(spike_idxs(i)),3),30,'r','filled');
end

end

