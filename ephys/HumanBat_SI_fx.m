function [SI_info] = HumanBat_SI_fx(batdate,logger,unit,flights_to_smush)
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

for bbb = 1:12
    unit=bbb;
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));

load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));

%% 1. Calculate 3d rate maps for flying
bin_size = 20; % cm to section room into

% a. Get room boundaries and bin into 10cmx10cmx10cm voxels
figure(); hold on; title("Room Boundaries, All Flights, No Rest"); for i=1:length(ciholas_flight_struct_resort); scatter3(ciholas_flight_struct_resort{i}.pos(:,1),ciholas_flight_struct_resort{i}.pos(:,2),ciholas_flight_struct_resort{i}.pos(:,3),3,'b','filled'); end; hold off;
room_bounds = [-290,290; -260,260; 1,230];
room_bounds_new_origin = [room_bounds(1,:)-room_bounds(1,1);room_bounds(2,:)-room_bounds(2,1);room_bounds(3,:)];
room_2d = zeros(round(room_bounds_new_origin(1:2,2)/bin_size)');

% b. Smush all the flights together (take out the resting time) 
if strcmp(flights_to_smush,'all')
    clip = 0;
    ciholas_flying_r = []; ephys_flying_r = [];
    for i=1:size(ciholas_flight_struct_resort,2)
        st = ciholas_flight_struct_resort{i}.fstart_trimmed+clip; ed = ciholas_flight_struct_resort{i}.fend_trimmed+clip;
        ciholas_flying_r = [ciholas_flying_r;ciholas_r(st:ed,:)];
    
        % Make binary vector of ephys for every flight
        e_bin = zeros(1,ciholas_flight_struct_resort{i}.ephys_trimmed_t);
        if ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
            e_bin(ciholas_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
        end
        ephys_flying_r = [ephys_flying_r,e_bin];
    end
elseif strcmp(flights_to_smush,'clus')
    clip = 0;
    ciholas_flying_r = []; ephys_flying_r = [];
    for i=1:size(ciholas_flight_struct_resort,2)
        if ciholas_flight_struct_resort{i}.fclus~=1
            st = ciholas_flight_struct_resort{i}.fstart_trimmed+clip; ed = ciholas_flight_struct_resort{i}.fend_trimmed+clip;
            ciholas_flying_r = [ciholas_flying_r;ciholas_r(st:ed,:)];
        
            % Make binary vector of ephys for every flight
            e_bin = zeros(1,ciholas_flight_struct_resort{i}.ephys_trimmed_t);
            if ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
                e_bin(ciholas_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
            end
            ephys_flying_r = [ephys_flying_r,e_bin];
        end
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
        ciholas_flying_r = []; ephys_flying_r = [];
        for i=1:size(ciholas_flight_struct_resort,2)
            if ciholas_flight_struct_resort{i}.tripod_landing == index_
                st = ciholas_flight_struct_resort{i}.fstart_trimmed+clip; ed = ciholas_flight_struct_resort{i}.fend_trimmed+clip;
                ciholas_flying_r = [ciholas_flying_r;ciholas_r(st:ed,:)];
            
                % Make binary vector of ephys for every flight
                e_bin = zeros(1,ciholas_flight_struct_resort{i}.ephys_trimmed_t);
                if ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
                    e_bin(ciholas_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
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
        ciholas_flying_r = []; ephys_flying_r = [];
        for i=1:size(ciholas_flight_struct_resort,2)
            if ciholas_flight_struct_resort{i}.tripod_takeoff == index_
                st = ciholas_flight_struct_resort{i}.fstart_trimmed+clip; ed = ciholas_flight_struct_resort{i}.fend_trimmed+clip;
                ciholas_flying_r = [ciholas_flying_r;ciholas_r(st:ed,:)];
            
                % Make binary vector of ephys for every flight
                e_bin = zeros(1,ciholas_flight_struct_resort{i}.ephys_trimmed_t);
                if ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
                    e_bin(ciholas_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
                end
                ephys_flying_r = [ephys_flying_r,e_bin];
            end
        end
    end
end

% Plot concatenated flights to check
figure(); hold on; title(strcat(flights_to_smush," flights")); 
scatter3(ciholas_r(:,1),ciholas_r(:,2),ciholas_r(:,3),0.5,[0.8 0.8 0.8]);
scatter3(ciholas_flying_r(:,1),ciholas_flying_r(:,2),ciholas_flying_r(:,3),6,'b','filled'); hold off;

% c. Get time-spent in each voxel (occupancy map)
ciholas_flying_r_rescale = [ciholas_flying_r(:,1)-room_bounds(1,1)*10, ciholas_flying_r(:,2)-room_bounds(2,1)*10];

% Add fake coordinates
ciholas_flying_r_rescale = [ciholas_flying_r_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
figure('Name','Time spent in each 20x20 voxel'); hist3(ciholas_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
[occMap,Xedges,Yedges,binX,binY] = histcounts2(ciholas_flying_r_rescale(:,1),ciholas_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
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
            temp_count = temp_count + sum(ephys_flying_r((i-1)*floor(length(ephys_flying_r)/length(ciholas_flying_r_rescale))+1:(i)*floor(length(ephys_flying_r)/length(ciholas_flying_r_rescale))));
        end
    end
    ephys_binCount(j) = temp_count;
end

%save(strcat(exp_data_path,'ephys/logger13/extracted_data/unit',num2str(unit),'_binCount_flights',flights_to_smush),'ephys_binCount','coord_pairs','coord_pairs_unique','occMap','binX','binY','-v7.3')

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
spike_map_raw(1,1)=0;
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

% % Plot heat map with cluster X flights overlaid (RELIES ON Ciholas_flights_Ephys_Analysis script output!)
% figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap with flightpaths ",num2str(clus))); hold on; imagesc(place_map_smoothed_with_NaN'); axis off; 
% scatter(ciholas_flying_r_rescale(:,1)/193,ciholas_flying_r_rescale(:,2)/193,0.25,'g');
% for k=1:length(Flight_Group_Matrix.Clusters{clus})
%     flightnum = og_idx{Flight_Group_Matrix.Clusters{clus}(k)};
%     scatter((ciholas_flight_struct_resort{flightnum}.pos(:,1)-room_bounds(1,1))/bin_size*2+1,(ciholas_flight_struct_resort{flightnum}.pos(:,2)-room_bounds(2,1))/bin_size*2-6,1,'r');
% end
% xlim([0 500]); ylim([0 475]);
% hold off;
% 
% % Plot the clusters to check
% figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap with flightpaths ",num2str(clus))); hold on; imagesc(spike_map_smoothed_with_occMap','Interpolation','bilinear');
% scatter(ciholas_flying_r_rescale(:,1)/bin_size*2+1,ciholas_flying_r_rescale(:,2)/bin_size*2-6,0.5,'w');
% for k=1:length(Kx_flights)
%     flightnum = og_idx{Kx_flights(k)};
%     scatter((ciholas_flight_struct_resort{flightnum}.pos(:,1)-room_bounds(1,1))/bin_size*2+1,(ciholas_flight_struct_resort{flightnum}.pos(:,2)-room_bounds(2,1))/bin_size*2-6,2,'r');
% end
% for k=1:length(Mx_flights)
%     flightnum = og_idx{Mx_flights(k)};
%     scatter((ciholas_flight_struct_resort{flightnum}.pos(:,1)-room_bounds(1,1))/bin_size*2+1,(ciholas_flight_struct_resort{flightnum}.pos(:,2)-room_bounds(2,1))/bin_size*2-6,2,'g');
% end
% xlim([0 500]); ylim([0 475]);
% title(strcat("All Cluster ",num2str(clus)," Flights to K (red) Flights to M (green)"));
% hold off;  


%% ======= Compute the following Place-Field Parameters: ===========
% (i) Sparsity, (ii) Coherence - currently un-used, (iii) Information_per_spike, (iv) Information_per_second, (v) Area of place-field, (vi) peak firing rate,
% is this a place cell (based on a criterion of information per spike > 0.4) : ====


% % -----  SPARSITY  (computed for the SMOOTHED field): -----
% % Sparsity = <r_i>^2 / <r_i^2> = sum( p_i * r_i )^2 / sum( p_i * r_i^2 )
% %    Where:
% %       r_i = firing rate in bin i ;
% %       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins.
% % See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).
% 
% r_i = place_map_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED with NaN Place Field
% p_i = occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
%     sum( occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
% r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
% sparsity = sum( p_i .* r_i )^2 / sum( p_i .* ( r_i .^2 ) ) ; % Sparsity
% 
% % -----  COHERENCE  (computed for the UN-SMOOTHED field!!!!!!!): -----
% % JANKY
% % Coherence = "first-order autocorrelation of the Place field =
% % Correlation (correlation coefficient) between the vector of firing rates, r_i, and the
% % average firing rates of the 8 nearest neighbors of bin i.  Formally:
% %    Coherence = corrcoef( r_i, mean(8 neighbours of r_i) )
% % See: Muller RU, Kubie JL, J. Neurosci 9, 4101-4110 (1989).
% 
% r_i_8neighbors_avg = place_map_raw ;  % Initialize
% for ii_y_bin = 2 : length(unique(binX)) - 1 , % Loop over x-bins, NOT INCL. EDGES
%     for ii_x_bin = 2 : length(unique(binY)) - 1  % Loop over y-bins, NOT INCL. EDGES
%         matrix_3x3_of_neighbors = ...
%             place_map_raw( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
%         num_neighbors_notNaN = sum( ~isnan( matrix_3x3_of_neighbors([1:4 6:9]) ) ); % Here I did NOT count the central bin itself
%         sum_including_the_central_bin = nansum( matrix_3x3_of_neighbors(:) );
%         if ( num_neighbors_notNaN > 0 ), % If there ARE non-NaN neighbors -- only then do the average
%             r_i_8neighbors_avg( ii_y_bin, ii_x_bin ) = ... % Remove the central bin -- and then average
%                 ( sum_including_the_central_bin - place_map_raw( ii_y_bin, ii_x_bin ) ) / ...
%                 num_neighbors_notNaN ;
%         end
%     end
% end
% r_i_8neighbors_avg = r_i_8neighbors_avg( idx_notNaN_PlaceField_un_smoothed_rate_map ); % Take only bins where the Original (unsmoothed) field was not-NaN
% r_i_UNsmoothed = place_map_raw( idx_notNaN_PlaceField_un_smoothed_rate_map ); % Here, Use the UN-SMOOTHED Place Field
% r_i_8neighbors_avg = r_i_8neighbors_avg(:) ; r_i_UNsmoothed = r_i_UNsmoothed(:) ; % Turn these into Column vectors
% r_i_8neighbors_avg(isnan(r_i_8neighbors_avg))=0; r_i_8neighbors_avg(isinf(r_i_8neighbors_avg))=0; r_i_UNsmoothed(isinf(r_i_UNsmoothed))=0;
% ccc = corrcoef( r_i_8neighbors_avg, r_i_UNsmoothed); % Correlation coefficient
% rrr = ccc(2,1);
% coherence_Ztransformed = 0.5 * log( (1+rrr)/(1-rrr) ); % Fisher's z-transform of r-values

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
f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap")); hold on; colorbar; imagesc(place_map_smoothed_with_NaN'); axis off; title(strcat("Unit ",num2str(unit)," 2D ratemap smooth with NaN. SI: ",num2str(information_per_spike))); hold off;
if ~exist(strcat(exp_data_path,"ephys/logger13/extracted_data/unit_",num2str(unit),"_shuffles"))
    mkdir(strcat(exp_data_path,"ephys/logger13/extracted_data/unit_",num2str(unit),"_shuffles"));
    saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap_',flights_to_smush),'jpg');
else
    saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap_',flights_to_smush),'jpg');
end

%% Red dots figure 
ff=figure(); hold on; title("Flight with ephys plotted on top");
for i=1:length(ciholas_flight_struct_resort)
    if strcmp(flights_to_smush,'clus')
        if ciholas_flight_struct_resort{i}.fclus==1
            continue
        else
            plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,1),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,2),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
                    xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
            scatter3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed,1),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed,2),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed,3),'g');
            % Transform ephys timestamps to ciholas time
            % Make artificial vector the time duration of ciholas
            flight_dur = (ciholas_flight_struct_resort{i}.fend_trimmed -  ciholas_flight_struct_resort{i}.fstart_trimmed);
            flight_vec = [1:flight_dur]./120;
            if ciholas_flight_struct_resort{i}.ephys_trimmed{unit} == 0
                disp(strcat("No unit firing on flight ",num2str(i)));
                continue;
            else
                for j=1:length(ciholas_flight_struct_resort{i}.ephys_trimmed{unit})
                    nearest_ephys_point = dsearchn(flight_vec',ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(j)/1e6);
                    scatter3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed+nearest_ephys_point,1),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed+nearest_ephys_point,2),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed+nearest_ephys_point,3),'r','filled');
                end
            end
        end
    elseif strcmp(flights_to_smush,'all')
        plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,1),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,2),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
                    xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        scatter3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed,1),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed,2),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed,3),'g');
        % Transform ephys timestamps to ciholas time
        % Make artificial vector the time duration of ciholas
        flight_dur = (ciholas_flight_struct_resort{i}.fend_trimmed -  ciholas_flight_struct_resort{i}.fstart_trimmed);
        flight_vec = [1:flight_dur]./120;
        if ciholas_flight_struct_resort{i}.ephys_trimmed{unit} == 0
            disp(strcat("No unit firing on flight ",num2str(i)));
            continue;
        else
            for j=1:length(ciholas_flight_struct_resort{i}.ephys_trimmed{unit})
                nearest_ephys_point = dsearchn(flight_vec',ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(j)/1e6);
                scatter3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed+nearest_ephys_point,1),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed+nearest_ephys_point,2),ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed+nearest_ephys_point,3),'r','filled');
            end
        end
    end
end
hold off;
saveas(ff,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap_reddots_',flights_to_smush),'jpg');

%% 5. Check if the shuffles already exist for that cell
% If they exist, run a t test to assess significance of that unit's SI
if exist(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI_',flights_to_smush,'.mat'))
    SI_info.shuffled = 1;
    shuffles = load(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI_',flights_to_smush,'.mat'));
    [h,p] = ttest2(information_per_spike,shuffles.information_per_spike_shuff);
    SI_info.h = h;
    SI_info.p = p;
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
    disp("Running shuffling!");
    num_shuffles = 500;

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
                    temp_count = temp_count + sum(ephys_flying_r_shift((i-1)*round(length(ephys_flying_r_shift)/length(ciholas_flying_r_rescale))+1:(i)*round(length(ephys_flying_r_shift)/length(ciholas_flying_r_rescale))));
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
        spike_map_raw(1,1)=0;
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
        
        place_map_smoothed_with_NaN = (place_map_smoothed.* idx_timespent_density)./idx_timespent_density;
        occMap_ts_thresh_smoothed_with_NaN = (occMap_ts_thresh_smoothed.* idx_timespent_density)./idx_timespent_density;
        
        idx_notNaN_PlaceField = find( ~isnan( place_map_smoothed_with_NaN ) ); % Find the indexes of non-NaN bins
        idx_isNaN_PlaceField = find( isnan( place_map_smoothed_with_NaN  ) ); % Find the indexes of NaN bins
        
        idx_notNaN_PlaceField_un_smoothed_rate_map = find( ~isnan( place_map_raw ) ); % Find the indexes of non-NaN bins

        % Compute SI for this shuffle
        r_i = place_map_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED Place Field
        p_i = occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
            sum( occMap_ts_thresh_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
        r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
        % % %             r = mean( r_i ) ;
        r = sum( r_i .* p_i );
        information_per_spike_shuff(lb) = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)
    
        f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap shuffle ",num2str(lb)," I: ",num2str(information_per_spike_shuff(lb))),'visible','off'); hold on; imagesc(place_map_smoothed_with_NaN'); title(strcat("Unit ",num2str(unit)," 2D ratemap shuffle ",num2str(lb)," I: ",num2str(information_per_spike_shuff(lb)))); hold off;
        
        if ~exist(strcat(exp_data_path,"ephys/logger13/extracted_data/unit_",num2str(unit),"_shuffles"))
            mkdir(strcat(exp_data_path,"ephys/logger13/extracted_data/unit_",num2str(unit),"_shuffles"));
            %saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap_shuffle_',num2str(lb)),'jpg');
        else
            %saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap_shuffle_',num2str(lb)),'jpg');
        end
    end
    save(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI_',flights_to_smush,'.mat'),'information_per_spike_shuff','-v7.3')

    disp("Doing Sig Test with new shuffled data!")
    shuffles = load(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI_',flights_to_smush,'.mat'));
    [h,p] = ttest2(information_per_spike,shuffles.information_per_spike_shuff);
    SI_info.h = h;
    SI_info.p = p;
    if h==1 & information_per_spike < mean(shuffles.information_per_spike_shuff) 
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
        SI_info.SI_category = 'SI < shuffle SI';
        disp("Unit has less spatial information than the mean of the shuffles.");
        ww = figure();  hold on; title(strcat("SI of shuffles and Unit ",num2str(unit)," (red dot) ",num2str(information_per_spike),". P = ",num2str(p))); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
       saveas(ww,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_sig_boxplot.jpg'));
    elseif h==1 & mean(shuffles.information_per_spike_shuff) < information_per_spike
        SI_info.SI_category = 'Spatially Tuned. SI > shuffle SI';
        SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
        disp("This Unit has significant spatial information.");
        disp(strcat("SI of unit ",num2str(unit),":",num2str(information_per_spike)));
        ww = figure();  hold on; title(strcat("SI of shuffles and Unit ",num2str(unit)," (red dot) ",num2str(information_per_spike),". P = ",num2str(p))); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
        saveas(ww,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_sig_boxplot.jpg'));
    else
       SI_info.mean_shuffle_information_per_spike = mean(shuffles.information_per_spike_shuff);
       SI_info.SI_category = 'None';
       disp("Not spatially tuned unit.");
       ww = figure();  hold on; title(strcat("SI of shuffles and Unit ",num2str(unit)," (red dot) ",num2str(information_per_spike),". P = ",num2str(p))); boxplot(shuffles.information_per_spike_shuff); scatter(1,information_per_spike,'r','filled'); ylim([0 max(shuffles.information_per_spike_shuff)+information_per_spike-max(shuffles.information_per_spike_shuff)+0.1]);
       saveas(ww,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_sig_boxplot.jpg'));
    end
end
end

end
