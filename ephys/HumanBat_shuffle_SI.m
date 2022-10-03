function [] = HumanBat_shuffle_SI(batdate,logger)

% Mass shuffle code
num_shuffles = 500; %batdate=220406;  logger=13;
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');

% Load in the ciholas data struct
if exist(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'))
    load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
else
    disp("Pruned, resorted, and trimmed struct not made!")
    return
end

% Load in the raw ciholas data
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));

% For every neuron, make the shuffles
for Neurons = 1:length(ciholas_flight_struct_resort{1}.ephys_trimmed)
    unit=Neurons;
    if ~exist(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles'))
        mkdir(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles'))
    end

    bin_size = 20; % cm to section room into

    % a. Get room boundaries and bin into 10cmx10cmx10cm voxels
    %figure(); hold on; title("Room Boundaries, All Flights, No Rest"); for i=1:100; plot3(ciholas_flight_struct_resort{i}.pos(:,1),ciholas_flight_struct_resort{i}.pos(:,2),ciholas_flight_struct_resort{i}.pos(:,3)); end; hold off;
    room_bounds = [-290,290; -260,260; 1,230];
    room_bounds_new_origin = [room_bounds(1,:)-room_bounds(1,1);room_bounds(2,:)-room_bounds(2,1);room_bounds(3,:)];
    room_2d = zeros(round(room_bounds_new_origin(1:2,2)/bin_size)');

    % b. Smush all the flights together (take out the resting time)
    clip = 0;
    ciholas_flying_r = []; ephys_flying_r = [];
    for i=1:size(ciholas_flight_struct_resort,2)
        st = ciholas_flight_struct_resort{i}.fstart_trimmed+clip; ed = ciholas_flight_struct_resort{i}.fend_trimmed+clip;
        ciholas_flying_r = [ciholas_flying_r;ciholas_r(st:ed,:)];

        % Make binary vector of ephys for eery flight
        e_bin = zeros(1,size(ciholas_flight_struct_resort{i}.ephys_trimmed_t,2));
        if ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(1) ~= 0
            e_bin(ciholas_flight_struct_resort{i}.ephys_trimmed{unit}) = 1;
        end
        ephys_flying_r = [ephys_flying_r,e_bin];
    end
    figure(); scatter3(ciholas_flying_r(:,1),ciholas_flying_r(:,2),ciholas_flying_r(:,3),1);

    % c. Get time-spent in each voxel
    ciholas_flying_r_rescale = [ciholas_flying_r(:,1)-room_bounds(1,1)*10, ciholas_flying_r(:,2)-room_bounds(2,1)*10];
    figure('Name','Time spent in each 20x20 voxel'); hist3(ciholas_flying_r_rescale,'NBins',[size(room_2d,1),size(room_2d,2)]);
    [occMap,Xedges,Yedges,binX,binY] = histcounts2(ciholas_flying_r_rescale(:,1),ciholas_flying_r_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);
    %figure('Name','Samples occupancy map'); imagesc(rot90(occMap));
    occMap_ts = occMap/120;  %figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
    occMap_ts_sm = imgaussfilt(occMap_ts);   %figure('Name','TS occupancy map, smoothed'); imagesc(rot90(occMap_ts_sm));
    occMap_rollout = reshape(occMap_ts,[1,size(occMap,1)*size(occMap,2)]);

    % e. For every XY-bin, get the # of spikes in that bin (1xn vector of
    % spikes where n is the cumulative length of occupied bins)

    coord_pairs = [binX,binY]; coord_pairs_unique = unique(coord_pairs,'rows');
    
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
                    temp_count = temp_count + sum(ephys_flying_r_shift((i-1)*round(length(ephys_flying_r_shift)/length(ciholas_flying_r_rescale))+1:(i)*round(length(ephys_flying_r_shift)/length(ciholas_flying_r_rescale))));
                end
            end
            ephys_binCount_shuff(j) = temp_count;
        end
        %save(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_',num2str(lb),'_binCounts.mat'),'ephys_binCount_shuff','-v7.3')
    
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
                        %spike_map_norm_sm(i,j) = ephys_binCount_norm_sm(m);
                        %spike_map_norm(i,j) = ephys_binCount_norm(m);
                        spike_map_raw(ubinX(i),ubinY(j)) = ephys_binCount_shuff(m);
                    end
                end
            end
        end
       
        
        % Calculate the raw place map using the absolute number of spikes per bin divided by
        % the occupancy in seconds in each bin
        place_map_raw = spike_map_raw./occMap_ts_thresh;
        %figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap raw")); hold on; imagesc(place_map_raw');%,'Interpolation','bilinear'); xlim([0 50]); ylim([0 47]);
        
        % Create gaussian kernel: 
        sigma = 1.5; % Hafting I guess?
        hsize = 5*round(sigma)+1;
        gaussian_kernel = fspecial('gaussian',hsize,sigma);
        %figure(); surf(gaussian_kernel); % If you want to see how your kernel looks like
        
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
        information_per_spike(lb) = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)
    
        f = figure('Name',strcat("Unit ",num2str(unit)," 2D ratemap shuffle ",num2str(lb)," I: ",num2str(information_per_spike(lb))),'visible','off'); hold on; imagesc(place_map_smoothed_with_NaN'); title(strcat("Unit ",num2str(unit)," 2D ratemap shuffle ",num2str(lb)," I: ",num2str(information_per_spike(lb)))); hold off;
        saveas(f,strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/',num2str(batdate),'_14650_PlaceMap_shuffle_',num2str(lb)),'jpg');
    end
    save(strcat(exp_data_path,'ephys/logger13/extracted_data/unit_',num2str(unit),'_shuffles/shuffle_SI.mat'),'information_per_spike','-v7.3')
end

end