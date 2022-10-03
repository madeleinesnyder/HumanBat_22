% Make dataset for Single Unit Linear Regression 

clear all;
batdate=220407;  logger=13;
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
units = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
E_cat = []; K_cat = []; M_cat = []; B_cat = []; B_other_cat = []; table1 = {}; P_cat = []; PO_cat = []; D_cat = []; C_cat = [];

% Load in ciholas_flight_struct_resort and B_ephys_data and ciholas_r and human_r
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'));
load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));
cortex_flights.trajectoriesContinous = cortex_flights.trajectoriesContinous';

% Sort the ciholas_flight_struct_resort into chronological order
cv = []; seconds = 3;
for i=1:length(ciholas_flight_struct_resort)
    cv(i) = ciholas_flight_struct_resort{i}.fstart_trimmed;
end
[im,imk] = sort(cv);
for i=1:length(ciholas_flight_struct_resort)
    ciholas_flight_struct_chron{i} = ciholas_flight_struct_resort{imk(i)};
end

% Make occMap for the ciholas bat
bin_size=20;
room_bounds = [-290,290; -260,260; 1,230];
room_bounds_new_origin = [room_bounds(1,:)-room_bounds(1,1);room_bounds(2,:)-room_bounds(2,1);room_bounds(3,:)];
room_2d = zeros(round(room_bounds_new_origin(1:2,2)/bin_size)');

%% Make training dataset (2) (At least one flight from each cluster and each location, -3s and +3s)
for nn=1:length(units)
    %E_cat = []; K_cat = []; M_cat = []; B_cat = []; B_other_cat = []; table1 = {}; P_cat = []; PO_cat = []; D_cat = []; C_cat = [];

unit=units(nn);
% Find appropriate flights
dataset2_flights = [];
for i=1:12
    if i==1
        cluster_flights = [];
        for j=1:length(ciholas_flight_struct_chron)
           if ciholas_flight_struct_chron{j}.fclus == i
               cluster_flights = [cluster_flights,j];
           end
        end
        cluster_flight_indexes = randperm(length(cluster_flights),round(length(cluster_flights)/2));
    else
        cluster_flights = [];
        for j=1:length(ciholas_flight_struct_chron)
            if ciholas_flight_struct_chron{j}.fclus == i
                cluster_flights = [cluster_flights,j];
            end
        end
        cluster_flight_indexes = randperm(length(cluster_flights),round(length(cluster_flights)));
    end
    dataset2_flights = [dataset2_flights,cluster_flights(cluster_flight_indexes)];
end

for i=1:length(dataset2_flights)
    ciholas_flight_struct_dataset2{i} = ciholas_flight_struct_chron{dataset2_flights(i)};
end

% Make sure it's sorted.   
cv_f=[];
for i=1:length(ciholas_flight_struct_dataset2)
    cv_f = [cv_f,ciholas_flight_struct_dataset2{i}.fstart_trimmed];
end
[im,imk] = sort(cv_f);
for i=1:length(ciholas_flight_struct_dataset2)
    ciholas_flight_struct_dataset2_chron{i} = ciholas_flight_struct_dataset2{imk(i)};
end

% Plot subset of flights used for taining the model
figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
for i=1:length(ciholas_flight_struct_dataset2_chron)
    scatter3(ciholas_flight_struct_dataset2_chron{i}.pos(:,1),ciholas_flight_struct_dataset2_chron{i}.pos(:,2),ciholas_flight_struct_dataset2_chron{i}.pos(:,3));
end
hold off;

%% Make the data vectors for this set of flights

for i=1:length(ciholas_flight_struct_dataset2_chron)

    temp_fstart = ciholas_flight_struct_dataset2_chron{i}.fstart_trimmed-seconds*120;
    %temp_fstart_onesec = ciholas_flight_struct_dataset2_chron{i}.fstart_trimmed-120;
    temp_fend = ciholas_flight_struct_dataset2_chron{i}.fend_trimmed+seconds*120;
    
    % ============== SELF POSITION ================== %     
    temp_flight = ciholas_r(temp_fstart:temp_fend,:);
    temp_flight_rescale = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10];

    % Add fake coordinates and roll out the occmap
    temp_flight_rescale_buff = [temp_flight_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
    %figure('Name','Time spent in each 20x20 voxel'); hist3(temp_flight_rescale_buff,'NBins',[size(room_2d,1),size(room_2d,2)]);
    [occMap,Xedges,Yedges,binX,binY] = histcounts2(temp_flight_rescale_buff(:,1),temp_flight_rescale_buff(:,2),[size(room_2d,1),size(room_2d,2)]);
    occMap_ts = occMap/120;  %figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
    occMap_rollout = reshape(occMap_ts,[1,size(occMap_ts,1)*size(occMap_ts,2)]);

    % Get the bin of the occMap_rollout for each timepoint in the flight
    % (MAYBE FISHY IF THERE IS A PROBLEM LATER ON)
    %bin_rollout_space = [1:(size(room_2d,1)*size(room_2d,2)+size(room_2d,1)*size(room_2d,2))];
    bin_rollout = [];
    for j=1:length(binX)
        bin_rollout(j,:) = [binX(j),binY(j)];
    end
    temp_bin_rollout = [];
    for j=1:length(bin_rollout)
        if bin_rollout(j,1) <= bin_rollout(j,2)
            rollout_space_index = bin_rollout(j,1)*bin_rollout(j,2);
        else
            rollout_space_index = size(room_2d,1)*size(room_2d,2) + bin_rollout(j,1)*bin_rollout(j,2);
        end
        temp_bin_rollout(j) = rollout_space_index;
    end

    dataset2.B{i} = temp_bin_rollout;
    % Downsample
    dataset2.B{i} = downsample(temp_bin_rollout,120);

    % ============== CLUSTER ================== % 
    dataset2.Cluster{i} = ciholas_flight_struct_dataset2_chron{i}.fclus;

    % For this window of time, find which tripods M and K are at 
    
    % Ground truth of tripods
    tripods = [-2.4,0.56,1.1;
               -1.6,2.07,1.1;
               2.07,1.6,0.9;
               2.26,0.1,0.9;
               2.1,-1.6,1;
               0.36,-2.2,0.9;
               -0.511, -0.265, 0.610];

    tripods = tripods*1000;
    cutoff=600;

    % Get mean M and K positions
    %mean_M_pos = mean(human_r(temp_fstart:temp_fend,:,4));
    %mean_K_pos = mean(human_r(temp_fstart:temp_fend,:,3));
    %M_distance_from_tripods = pdist2(tripods,mean_M_pos, 'euclidean');
    %K_distance_from_tripods = pdist2(tripods,mean_K_pos, 'euclidean');
    %for j=1:length(M_distance_from_tripods)
    %dataset2.M_tripod{i} = find(M_distance_from_tripods == min(M_distance_from_tripods));
    %end
    %figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]);
    %scatter3(human_r(temp_fstart:temp_fend,1,4),human_r(temp_fstart:temp_fend,2,4),human_r(temp_fstart:temp_fend,3,4),'b');
    %scatter3(mean_M_pos(1,1),mean_M_pos(1,2),mean_M_pos(1,3),100,'black','filled'); hold off;
    
    %for j=1:length(K_distance_from_tripods)
    %dataset2.K_tripod{i} = find(K_distance_from_tripods == min(K_distance_from_tripods));
    %end
    %figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]);
    %scatter3(human_r(temp_fstart:temp_fend,1,3),human_r(temp_fstart:temp_fend,2,3),human_r(temp_fstart:temp_fend,3,3),'b');
    %scatter3(mean_K_pos(1,1),mean_K_pos(1,2),mean_K_pos(1,3),100,'black','filled'); hold off;
    %dataset2.Cluster{i} = ciholas_flight_struct_dataset2_chron{i}.fclus;

    % ============== EPHYS ================== % 
    flight_ephys_B = {}; ephys=0;   tt= {}; rr={};
    f_start_seconds = temp_fstart/120;  f_end_seconds = temp_fend/120;   % Get ephys samples in seconds
    % Find the closest ephys sample to f_start_seconds*1e6 
    ephys_f_start = f_start_seconds*1e6; % This is the start time for that ephys chunk. Cut out that ephys chunk
    ephys_f_end = f_end_seconds*1e6;

    % Make a vector tt that is the length of the seconds 1e6 in
    % that flight
    tt = [0:(ephys_f_end-ephys_f_start)];   
    temp_rr = zeros(length(tt),1);

    % Fill in the rest of the timestamps with ephys_f_start
    % subtracted for that ephys stretch
    [B,ephys_idxs_in_range] = find(ephys_f_start <= B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps <= ephys_f_end);
    dataset2.E_check{i} = ephys_idxs_in_range;
    E_temp = {}; E_temp_summed = [];
    for j=1:length(dataset2.B{i})
        [B,E_in_sec] = find(ephys_f_start+(j-1)*1e6 <= B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps <= ephys_f_start+(j)*1e6);
        if isempty(E_in_sec)
            E_temp{j} = 0;
        else
            E_temp{j} = E_in_sec;
        end
        for j=1:length(E_temp)
            E_temp_summed(j) = length(E_temp{j});
        end
        dataset2.E{i} = E_temp_summed;
    end
        %dataset2.E{i} = 0; % zeros(1,length(tt));
        %dataset2.E_samples{i} = length(tt);
    %else
    %    for j=1:length(ephys_idxs_in_range)
    %        temp_ts = B_ephys_data.TT_unit(unit).AlignedTimestamps(ephys_idxs_in_range(j));
     %       temp_rr(round(temp_ts-ephys_f_start)) = 1;
     %   end
        %dataset2.E{i} = find(temp_rr==1);
        %dataset2.E_samples{i} = length(tt);
   % end

    % ============== HUMANS ================== % 

    % KEVIN
    temp_flight = [];
    temp_flight = human_r(temp_fstart:temp_fend,:,3);
    temp_flight_rescale = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10];

    % Add fake coordinates and roll out the occmap
    temp_flight_rescale_buff = [temp_flight_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
    %figure('Name','Time spent in each 20x20 voxel'); hist3(temp_flight_rescale_buff,'NBins',[size(room_2d,1),size(room_2d,2)]);
    [occMap,Xedges,Yedges,binX,binY] = histcounts2(temp_flight_rescale_buff(:,1),temp_flight_rescale_buff(:,2),[size(room_2d,1),size(room_2d,2)]);
    occMap_ts = occMap/120;  %figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
    occMap_rollout = reshape(occMap_ts,[1,size(occMap_ts,1)*size(occMap_ts,2)]);

    % Get the bin of the occMap_rollout for each timepoint in the flight
    bin_rollout = [];
    for j=1:length(binX)
        bin_rollout(j,:) = [binX(j),binY(j)];
    end
    temp_bin_rollout = [];
    for j=1:length(bin_rollout)
        if bin_rollout(j,1) <= bin_rollout(j,2)
            rollout_space_index = bin_rollout(j,1)*bin_rollout(j,2);
        else
            rollout_space_index = size(room_2d,1)*size(room_2d,2) + bin_rollout(j,1)*bin_rollout(j,2);
        end
        temp_bin_rollout(j) = rollout_space_index;
    end

    dataset2.K{i} = temp_bin_rollout;
    dataset2.K{i} = downsample(temp_bin_rollout,120);

    % MADELEINE
    temp_flight = [];
    temp_flight = human_r(temp_fstart:temp_fend,:,4);
    temp_flight_rescale = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10];

    % Add fake coordinates and roll out the occmap
    temp_flight_rescale_buff = [temp_flight_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
    %figure('Name','Time spent in each 20x20 voxel'); hist3(temp_flight_rescale_buff,'NBins',[size(room_2d,1),size(room_2d,2)]);
    [occMap,Xedges,Yedges,binX,binY] = histcounts2(temp_flight_rescale_buff(:,1),temp_flight_rescale_buff(:,2),[size(room_2d,1),size(room_2d,2)]);
    occMap_ts = occMap/120;  %figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
    occMap_rollout = reshape(occMap_ts,[1,size(occMap_ts,1)*size(occMap_ts,2)]);

    % Get the bin of the occMap_rollout for each timepoint in the flight
    bin_rollout = [];
    for j=1:length(binX)
        bin_rollout(j,:) = [binX(j),binY(j)];
    end
    temp_bin_rollout = [];
    for j=1:length(bin_rollout)
        if bin_rollout(j,1) <= bin_rollout(j,2)
            rollout_space_index = bin_rollout(j,1)*bin_rollout(j,2);
        else
            rollout_space_index = size(room_2d,1)*size(room_2d,2) + bin_rollout(j,1)*bin_rollout(j,2);
        end
        temp_bin_rollout(j) = rollout_space_index;
    end

    dataset2.M{i} = temp_bin_rollout;
    dataset2.M{i} = downsample(temp_bin_rollout,120);

    % ============== PDPO ================== % 
    pre_ = zeros(length(dataset2.E{i}),1); p_i = [1,2,3]; pre_(p_i) = 1; 
    post_ = zeros(length(dataset2.E{i}),1); po_i = [length(dataset2.E{i})-2,length(dataset2.E{i})-1,length(dataset2.E{i})]; post_(po_i) = 1;
    dur_ = (~pre_ & ~post_) ;%zeros(length(dataset2.E{i}),1); d_i = [length(dataset2.E{i})-3-2,length(dataset2.E{i})-3-1,length(dataset2.E{i})-3]; dur_(d_i) = 1;

    dataset2.pre{i} = pre_; dataset2.dur{i} = dur_; dataset2.post{i} = post_;

    % ============== OTHER BAT ================== % 
    temp_flight = cortex_flights.trajectoriesContinous(temp_fstart:temp_fend,:)*1000;
    temp_flight_rescale = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10];

    % Add fake coordinates and roll out the occmap
    temp_flight_rescale_buff = [temp_flight_rescale;[0 0]; [room_bounds_new_origin(1,2)*10 room_bounds_new_origin(2,2)*10]];
    %figure('Name','Time spent in each 20x20 voxel'); hist3(temp_flight_rescale_buff,'NBins',[size(room_2d,1),size(room_2d,2)]);
    [occMap,Xedges,Yedges,binX,binY] = histcounts2(temp_flight_rescale_buff(:,1),temp_flight_rescale_buff(:,2),[size(room_2d,1),size(room_2d,2)]);
    occMap_ts = occMap/120;  %figure('Name','TS occupancy map'); imagesc(rot90(occMap_ts));
    occMap_rollout = reshape(occMap_ts,[1,size(occMap_ts,1)*size(occMap_ts,2)]);

    % Get the bin of the occMap_rollout for each timepoint in the flight
    bin_rollout_space = [1:(size(room_2d,1)*size(room_2d,2)+size(room_2d,1)*size(room_2d,2))];
    bin_rollout = []; temp_bin_rollout = [];
    for j=1:length(binX)
        bin_rollout(j,:) = [binX(j),binY(j)];
    end
    for j=1:length(bin_rollout)
        if bin_rollout(j,1) <= bin_rollout(j,2)
            rollout_space_index = bin_rollout(j,1)*bin_rollout(j,2);
        else
            rollout_space_index = size(room_2d,1)*size(room_2d,2) + bin_rollout(j,1)*bin_rollout(j,2);
        end
        temp_bin_rollout(j) = rollout_space_index;
    end

    dataset2.B_other{i} = temp_bin_rollout;  
    dataset2.B_other{i} = downsample(temp_bin_rollout,120);
end

%% Table 1:

%E_cat = []; K_cat = []; M_cat = []; B_cat = []; B_other_cat = [];
for i=1:length(dataset2.B)

    % Concatenate the ephys data
%     E_zeros = zeros(1,dataset2.E_samples{i});
%     try
%         E_zeros(dataset2.E{i}) = 1;
%     catch
%         disp("No spikes for this flight");
%     end
%    to_trim = mod(length(E_zeros),length(dataset2.B{i}));
%    E_zeros = E_zeros(to_trim+1:end)';
    E_cat = [E_cat; dataset2.E{i}'];

    % Concatenate and resample the position data (LAGGING DUE TO ROUNDING)
%     B_resampled = [];
%     resample_factor = length(E_zeros)/length(dataset2.B{i});
%     for j=1:length(dataset2.B{i})
%         B_resampled = [B_resampled;repmat(dataset2.B{i}(j),1,resample_factor)'];
%     end
    B_cat = [B_cat;dataset2.B{i}'];

%     B_other_resampled = [];
%     resample_factor_o = length(E_zeros)/length(dataset2.B_other{i});
%     for j=1:length(dataset2.B_other{i})
%         B_other_resampled = [B_other_resampled;repmat(dataset2.B_other{i}(j),1,resample_factor_o)'];
%     end
    B_other_cat = [B_other_cat;dataset2.B_other{i}'];

    K_cat = [K_cat;dataset2.K{i}'];
    M_cat = [M_cat;dataset2.M{i}'];

    % Concatenate and resample the M and K data 
    %K_resampled = []; M_resampled = [];
    %resample_factor_KM = length(E_zeros)/length(dataset2.M_tripod{i});
    %M_resampled = [M_resampled; repmat(dataset2.M_tripod{i},1,length(dataset2.B{i}))'];
    %K_resampled = [K_resampled; repmat(dataset2.K_tripod{i},1,length(dataset2.B{i}))'];
    %M_cat = [M_cat;M_resampled]; K_cat = [K_cat;K_resampled];

    % Label whether in pre, during, or post positions
    P_cat = [P_cat; dataset2.pre{i}];
    D_cat = [D_cat; dataset2.dur{i}];
    PO_cat = [PO_cat; dataset2.post{i}];

    clus_resampled = [];
    clus_resampled = [clus_resampled; repmat(dataset2.Cluster{i},1,length(dataset2.B{i}))'];
    C_cat = [C_cat; clus_resampled];
end

% Make a table
%ID_cat = ones(length(E_cat),1)*unit;
table_ = [E_cat,B_cat,B_other_cat,K_cat,M_cat,C_cat,P_cat,D_cat,PO_cat];
ID_cat = []; K_cat=[]; B_cat=[];M_cat=[]; K_resampled=[]; M_resampled=[]; B_other_resampled=[]; E_cat=[]; B_other_cat = []; C_cat = []; P_cat = []; D_cat = []; PO_cat = [];

save(strcat(exp_data_path,'GLM/training/',num2str(batdate),'_14650_dataset2_allFlights_',num2str(nn),'.mat'),'table_','-v7.3');
end
