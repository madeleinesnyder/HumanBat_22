% Make dataset for Single Unit Linear Regression with 100ms bins and rv
% integer Fr. Include the distance of the bat to each human as variables
% and the positio of self other bat and humans in coordinates based on the
% lower left corner of the room. WIth -1s and +1s window

clear all;
only_pre = 0;
only_dur = 0;

batdate=220408;  logger=15;
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
E_cat = []; P_cat = []; PO_cat = []; D_cat = []; C_cat = [];
Kx_cat = []; Mx_cat = []; Bx_cat = []; Bx_other_cat = []; 
Ky_cat = []; My_cat = []; By_cat = []; By_other_cat = []; 
BMX_cat = []; BMY_cat = []; BKX_cat = []; BKY_cat = []; BBOX_cat = []; BBOY_cat = [];
table1 = {};

% Load in ciholas_flight_struct_resort and B_ephys_data and ciholas_r and human_r
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
units = [1:length(B_ephys_data.TT_unit)];
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'));
load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));
cortex_flights.trajectoriesContinous = cortex_flights.trajectoriesContinous';

% Sort the ciholas_flight_struct_resort into chronological order
cv = []; seconds = 1;
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
E_cat = []; P_cat = []; PO_cat = []; D_cat = []; C_cat = [];
Kx_cat = []; Mx_cat = []; Bx_cat = []; Bx_other_cat = []; 
Ky_cat = []; My_cat = []; By_cat = []; By_other_cat = []; 
BMX_cat = []; BMY_cat = []; BKX_cat = []; BKY_cat = []; BBOX_cat = []; BBOY_cat = [];
table1 = {};

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
    temp_fend = ciholas_flight_struct_dataset2_chron{i}.fend_trimmed+seconds*120;
    
    % ============== SELF POSITION X Y ================== %     
    temp_flight = ciholas_r(temp_fstart:temp_fend,:);
    temp_flight_rescale = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10];
    temp_flight_rescale_AD = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10,temp_flight(:,3)*10];
    [occMap,Xedges,Yedges,binX,binY] = histcounts2(temp_flight_rescale(:,1),temp_flight_rescale(:,2),[size(room_2d,1),size(room_2d,2)]);

    % Downsample flight data to 100ms bins 
    dataset2.B_x{i} = downsample(binX,12);%downsample(temp_flight_rescale(:,1),12);
    dataset2.B_y{i} = downsample(binY,12);%downsample(temp_flight_rescale(:,2),12);

    % ============== CLUSTER ================== % 
    dataset2.Cluster{i} = ciholas_flight_struct_dataset2_chron{i}.fclus;

    
    % ============== HUMANS ================== % 

    % KEVIN
    temp_flightK = [];
    temp_flightK = human_r(temp_fstart:temp_fend,:,3);
    temp_flight_rescaleK = [temp_flightK(:,1)-room_bounds(1,1)*10, temp_flightK(:,2)-room_bounds(2,1)*10];
    [occMap,Xedges,Yedges,binXK,binYK] = histcounts2(temp_flight_rescaleK(:,1),temp_flight_rescaleK(:,2),[size(room_2d,1),size(room_2d,2)]);

    % Downsample flight data to 100ms bins 
    dataset2.K_x{i} = downsample(binXK,12);%downsample(temp_flight_rescaleK(:,1),12);
    dataset2.K_y{i} = downsample(binYK,12);%downsample(temp_flight_rescaleK(:,2),12);

   
    % MADELEINE
    temp_flightM = [];
    temp_flightM = human_r(temp_fstart:temp_fend,:,4);
    temp_flight_rescaleM = [temp_flightM(:,1)-room_bounds(1,1)*10, temp_flightM(:,2)-room_bounds(2,1)*10];
    [occMap,Xedges,Yedges,binXM,binYM] = histcounts2(temp_flight_rescaleM(:,1),temp_flight_rescaleM(:,2),[size(room_2d,1),size(room_2d,2)]);

    % Downsample flight data to 100ms bins 
    dataset2.M_x{i} = downsample(binXM,12);%downsample(temp_flight_rescaleM(:,1),12);
    dataset2.M_y{i} = downsample(binYM,12);%downsample(temp_flight_rescaleM(:,2),12);



    % ============== OTHER BAT ================== % 
    
    temp_flight = cortex_flights.trajectoriesContinous(temp_fstart:temp_fend,:)*1000;
    temp_flight_rescaleBO = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10];
    temp_flight_rescaleBO_AD = [temp_flight(:,1)-room_bounds(1,1)*10, temp_flight(:,2)-room_bounds(2,1)*10,temp_flight(:,3)*10];
    [occMap,Xedges,Yedges,binXBO,binYBO] = histcounts2(temp_flight_rescaleBO(:,1),temp_flight_rescaleBO(:,2),[size(room_2d,1),size(room_2d,2)]);

    % Downsample flight data to 100ms bins 
    dataset2.Bx_other{i} = downsample(binXBO,12);%downsample(temp_flight_rescale_O(:,1),12);
    dataset2.By_other{i} = downsample(binYBO,12);%downsample(temp_flight_rescale_O(:,2),12);

   


    % ============== BAT'S DISTANCES FROM M AND K and Other Bat ================== % 

    M_pos = human_r(temp_fstart:temp_fend,:,4);     K_pos = human_r(temp_fstart:temp_fend,:,3);     
    B_pos = ciholas_r(temp_fstart:temp_fend,:);
    B_pos_rs = [B_pos(:,1)-room_bounds(1,1)*10, B_pos(:,2)-room_bounds(2,1)*10];
    M_pos_rs = [M_pos(:,1)-room_bounds(1,1)*10, M_pos(:,2)-room_bounds(2,1)*10];       K_pos_rs = [K_pos(:,1)-room_bounds(1,1)*10, K_pos(:,2)-room_bounds(2,1)*10]; 
    B2M = [];  B2K = [];  B2OB = [];
    for b = 1:length(B_pos_rs)
        B2M(b,:) = pdist2(B_pos_rs(b,:),M_pos_rs(b,:),'euclidean');
        B2K(b,:) = pdist2(B_pos_rs(b,:),K_pos_rs(b,:),'euclidean');
        B2OB(b,:) = pdist2(B_pos_rs(b,:),temp_flight_rescaleBO(b,:),'euclidean');
    end
    dataset2.B2M{i} = downsample(B2M,12);
    dataset2.B2K{i} = downsample(B2K,12);
    dataset2.B2OB{i} = downsample(B2OB,12);


    % ============== BAT'S DISTANCES FROM Home and Goal ================== % 

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
    Bat_distance_from_tripods = pdist2(tripods,mean(temp_flight_rescale_AD(end-5:end,:)), 'euclidean');
    goal_tripod = find(Bat_distance_from_tripods == min(Bat_distance_from_tripods));

    % G_factor is 1 if rewarded flight, and 2 if rest flight 
    if goal_tripod >= 8
        Home = mean(temp_flight_rescale(end-5:end,:));
        Goal =  mean(temp_flight_rescale(end-5:end,:));
        dataset2.Home{i} = mean(temp_flight_rescale(end-5:end,:));
        dataset2.Goal{i} = mean(temp_flight_rescale(end-5:end,:));
        dataset2.D2Home{i} = pdist2(temp_flight_rescale,Home,'euclidean');
        dataset2.D2Goal{i} = pdist2(temp_flight_rescale,Goal,'euclidean');
        dataset2.G_factor{i} = 2;
    else
        Goal =  mean(temp_flight_rescale(end-5:end,:));
        dataset2.Goal{i} = mean(temp_flight_rescale(end-5:end,:));
        dataset2.D2Home{i} = pdist2(temp_flight_rescale,Home,'euclidean')';
        dataset2.D2Goal{i} = pdist2(temp_flight_rescale,Goal,'euclidean');
        dataset2.G_factor{i} = 1;
    end

    % Downsample
    dataset2.D2Home{i} = downsample(dataset2.D2Home{i},12);
    dataset2.D2Goal{i} = downsample(dataset2.D2Goal{i},12);


    % ============== Factor that is whether flight is going to M or K Joint Rest or Solo Rest ================== % 

    M_pos = human_r(temp_fstart:temp_fend,:,4);     K_pos = human_r(temp_fstart:temp_fend,:,3);     OB_pos = temp_flight_rescaleBO_AD;
    cutoff = 600;
    %BD_M = pdist2(mean(temp_flight_rescale_AD(end-3:end,:)),mean(M_pos(end-3:end,:)));
    %BD_K = pdist2(mean(temp_flight_rescale_AD(end-3:end,:)),mean(K_pos(end-3:end,:)));
    %BD_OB = pdist2(mean(temp_flight_rescale_AD(end-3:end,:)),mean(OB_pos(end-3:end,:)));
    BD_M = pdist2(mean(ciholas_flight_struct_dataset2_chron{i}.pos(end-3:end,:)),mean(ciholas_flight_struct_dataset2_chron{i}.pos_M(end-3:end,:)));
    BD_K = pdist2(mean(ciholas_flight_struct_dataset2_chron{i}.pos(end-3:end,:)),mean(ciholas_flight_struct_dataset2_chron{i}.pos_K(end-3:end,:)));
    BD_OB = pdist2(mean(ciholas_flight_struct_dataset2_chron{i}.pos(end-3:end,:)),mean(ciholas_flight_struct_dataset2_chron{i}.pos_cortex_bat(:,end-3:end)'));
    BD_vec = [BD_M, BD_K, BD_OB];

    % Dataset.Who is 1 if Madeliene, 2 if Kevin, 3 if other bat, 4 if alone
    if any(BD_vec < cutoff)
        dataset2.Who_Landing{i} = find(BD_vec == min(BD_vec));
    else
        dataset2.Who_Landing{i} = 4;
    end
    


    % ============== Factor that is whether flight is coming from M or K ================== % 
    BD_vec = [];
    M_pos = human_r(temp_fstart:temp_fend,:,4);     K_pos = human_r(temp_fstart:temp_fend,:,3);     OB_pos = temp_flight_rescaleBO_AD;
    cutoff = 600;
    %BD_M = pdist2(mean(temp_flight_rescale_AD(end-3:end,:)),mean(M_pos(end-3:end,:)));
    %BD_K = pdist2(mean(temp_flight_rescale_AD(end-3:end,:)),mean(K_pos(end-3:end,:)));
    %BD_OB = pdist2(mean(temp_flight_rescale_AD(end-3:end,:)),mean(OB_pos(end-3:end,:)));
    %BD_vec = [BD_M, BD_K, BD_OB];
    BD_M = pdist2(mean(ciholas_flight_struct_dataset2_chron{i}.pos(end-3:end,:)),mean(ciholas_flight_struct_dataset2_chron{i}.pos_M(end-3:end,:)));
    BD_K = pdist2(mean(ciholas_flight_struct_dataset2_chron{i}.pos(end-3:end,:)),mean(ciholas_flight_struct_dataset2_chron{i}.pos_K(end-3:end,:)));
    BD_OB = pdist2(mean(ciholas_flight_struct_dataset2_chron{i}.pos(end-3:end,:)),mean(ciholas_flight_struct_dataset2_chron{i}.pos_cortex_bat(:,end-3:end)'));
    BD_vec = [BD_M, BD_K, BD_OB];

    % Dataset.Who is 1 if Madeliene, 2 if Kevin, 3 if other bat, 4 if alone
    if any(BD_vec < cutoff)
        dataset2.Who_Takeoff{i} = find(BD_vec == min(BD_vec));
    else
        dataset2.Who_Takeoff{i} = 4;
    end
    

    % Get mean M and K positions
%     M_pos = (human_r(temp_fstart:temp_fend,:,4));
%     Mposx = squeeze(temp_flight_rescaleM(:,1)); Bposx = squeeze(temp_flight_rescale(:,1));      Mposy = squeeze(temp_flight_rescaleM(:,2)); Bposy = squeeze(temp_flight_rescale(:,2));
%     BMX  = abs(Mposx - Bposx);   BMY  = abs(Mposy - Bposy);
%     dataset2.BMX{i} = downsample(BMX,12);
%     dataset2.BMY{i} = downsample(BMY,12);
% 
%     K_pos = (human_r(temp_fstart:temp_fend,:,3));
%     Kposx = squeeze(temp_flight_rescaleK(:,1));       Kposy = squeeze(temp_flight_rescaleK(:,2)); 
%     BKX  = abs(Kposx - Bposx);   BKY  = abs(Kposy - Bposy);
%     dataset2.BKX{i} = downsample(BKX,12);
%     dataset2.BKY{i} = downsample(BKY,12);
% 
%     BOposx = squeeze(temp_flight_rescaleBO(:,1));       BOposy = squeeze(temp_flight_rescaleBO(:,2)); 
%     BBOX  = abs(BOposx - Bposx);   BBOY  = abs(BOposy - Bposy);
%     dataset2.BBOX{i} = downsample(BBOX,12);
%     dataset2.BBOY{i} = downsample(BBOY,12);



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

    [B,ephys_idxs_in_range] = find(ephys_f_start <= B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps <= ephys_f_end);
    dataset2.E_check{i} = ephys_idxs_in_range;
    E_temp = {}; E_temp_summed = [];
    for j=1:length(dataset2.B2K{i})
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




    % ============== SANITY CHECK; DISTANCE FROM PEAK FIRING LOCATION ================== % 

    [B,ephys_idxs_in_range] = find(ephys_f_start <= B_ephys_data.TT_unit(unit).AlignedTimestamps & B_ephys_data.TT_unit(unit).AlignedTimestamps <= ephys_f_end);
    Peak_FR_idxs = find(dataset2.E{i} == max(dataset2.E{i}));
    Peak_FR_idx = Peak_FR_idxs(1);
    ds_tf = downsample(temp_flight_rescale_AD,12);
    Peak_Loc = ds_tf(Peak_FR_idx,:);
    dataset2.DFPF{i} = pdist2(temp_flight_rescale,Peak_Loc(1:2),'euclidean');
    dataset2.DFPF{i} = downsample(dataset2.DFPF{i},12);

    % ============== PDPO ================== % 
    pre_ = zeros(length(dataset2.E{i}),1); p_i = [1:10]; pre_(p_i) = 1; 
    post_ = zeros(length(dataset2.E{i}),1); po_i = [length(dataset2.E{i})-9:length(dataset2.E{i})]; post_(po_i) = 3;
    dur_ = (~pre_ & ~post_)*2; %zeros(length(dataset2.E{i}),1); d_i = [length(dataset2.E{i})-3-2,length(dataset2.E{i})-3-1,length(dataset2.E{i})-3]; dur_(d_i) = 1;
    PDP = pre_+dur_+post_;

    dataset2.PDP{i} = PDP;



    % ============== IF ONLY PRE ================== % 
    if only_pre == 1
        dataset2.B2K{i} = dataset2.B2K{i}(1:11);
        dataset2.B2M{i} = dataset2.B2M{i}(1:11);
        dataset2.E{i} = dataset2.E{i}(1:11);
        dataset2.D2Home{i} = dataset2.D2Home{i}(1:11);
        dataset2.D2Goal{i} = dataset2.D2Goal{i}(1:11);
        dataset2.DFPF{i} = dataset2.DFPF{i}(1:11);
    elseif only_dur ==1
        dataset2.B2K{i} = dataset2.B2K{i}(dur_ == 2);
        dataset2.B2M{i} = dataset2.B2M{i}(dur_ == 2);
        dataset2.E{i} = dataset2.E{i}(dur_ == 2);
        dataset2.D2Home{i} = dataset2.D2Home{i}(dur_ == 2);
        dataset2.D2Goal{i} = dataset2.D2Goal{i}(dur_ == 2);
        dataset2.DFPF{i} = dataset2.DFPF{i}(dur_ == 2);
    end
    

end

%% Table 1:

E_cat = []; P_cat = []; PO_cat = []; D_cat = []; C_cat = [];
%Kx_cat = []; Mx_cat = []; Bx_cat = []; Bx_other_cat = []; 
%Ky_cat = []; My_cat = []; By_cat = []; By_other_cat = []; 
%BMX_cat = []; BMY_cat = []; BKX_cat = []; BKY_cat = []; BBOX_cat = []; BBOY_cat = [];
D2Goal_cat = []; D2Home_cat = []; Home_cat = []; Goal_cat = []; Who_Takeoff_cat = []; Who_Landing_cat = []; D2DFPF_cat = []; G_factor_cat = [];
B2K_cat = []; B2M_cat = []; PDP_cat = []; B2OB_cat = [];
table_ = [];
for i=1:length(dataset2.B2K)

    E_cat = [E_cat; dataset2.E{i}'];

    B2K_cat = [B2K_cat;dataset2.B2K{i}];
    B2M_cat = [B2M_cat;dataset2.B2M{i}];
    B2OB_cat = [B2OB_cat;dataset2.B2OB{i}];

    D2Goal_cat = [D2Goal_cat;dataset2.D2Goal{i}];
    D2Home_cat = [D2Home_cat;dataset2.D2Home{i}];
    D2DFPF_cat = [D2DFPF_cat;dataset2.DFPF{i}];


%     Bx_cat = [Bx_cat;dataset2.B_x{i}];
%     By_cat = [By_cat;dataset2.B_y{i}];
% 
%     Bx_other_cat = [Bx_other_cat;dataset2.Bx_other{i}];
%     By_other_cat = [By_other_cat;dataset2.By_other{i}];
% 
%     Kx_cat = [Kx_cat;dataset2.K_x{i}];
%     Mx_cat = [Mx_cat;dataset2.M_x{i}];
%     Ky_cat = [Ky_cat;dataset2.K_y{i}];
%     My_cat = [My_cat;dataset2.M_y{i}];

%     BMX_cat = [BMX_cat;dataset2.BMX{i}];
%     BMY_cat = [BMY_cat;dataset2.BMY{i}];
%     BKX_cat = [BKX_cat;dataset2.BKX{i}];
%     BKY_cat = [BKY_cat;dataset2.BKY{i}];
%     BBOX_cat = [BBOX_cat;dataset2.BBOX{i}];
%     BBOY_cat = [BBOY_cat;dataset2.BBOY{i}];

%     P_cat = [P_cat; dataset2.pre{i}];
%     D_cat = [D_cat; dataset2.dur{i}];
%     PO_cat = [PO_cat; dataset2.post{i}];
    if only_pre == 0
        PDP_cat = [PDP_cat;dataset2.PDP{i}];
    end
    
    C_cat = [C_cat; repmat(dataset2.Cluster{i},1,length(dataset2.B2K{i}))'];
    %Goal_cat = [Goal_cat; repmat(dataset2.Goal{i},1,length(dataset2.B2K{i}))'];
    %Home_cat = [Home_cat; repmat(dataset2.Home{i},1,length(dataset2.B2K{i}))'];
    G_factor_cat = [G_factor_cat; repmat(dataset2.G_factor{i},1,length(dataset2.B2K{i}))'];
    Who_Landing_cat = [Who_Landing_cat; repmat(dataset2.Who_Landing{i},1,length(dataset2.B2K{i}))'];
    Who_Takeoff_cat = [Who_Takeoff_cat; repmat(dataset2.Who_Takeoff{i},1,length(dataset2.B2K{i}))'];
end

% Make a table
if only_pre == 1 | only_dur == 1
    table_ = [E_cat,B2M_cat,B2K_cat,B2OB_cat,B2G_cat,B2H_cat,C_cat];
else
    table_ = [E_cat,B2M_cat,B2K_cat,B2OB_cat,D2Goal_cat,D2Home_cat, Who_Landing_cat, Who_Takeoff_cat, G_factor_cat, C_cat,D2DFPF_cat,PDP_cat];
end
%table_ = [E_cat,B2M_cat,B2K_cat,B2OB_cat,B2G_cat,B2H_cat,C_cat,PDP_cat];
%table_ = [E_cat,Bx_cat,By_cat,BMX_cat,BMY_cat,BKX_cat,BKY_cat,BBOX_cat,BBOY_cat,C_cat,P_cat,D_cat,PO_cat];
%ID_cat = []; K_cat=[]; B_cat=[];M_cat=[]; K_resampled=[]; M_resampled=[]; B_other_resampled=[]; E_cat=[]; B_other_cat = []; C_cat = []; P_cat = []; D_cat = []; PO_cat = []; DF_K = []; DF_M = [];

if only_pre == 1
    save(strcat(exp_data_path,'GLM/training/',num2str(batdate),'_14650_datasetDistanceOnlyPre_allFlights_',num2str(nn),'.mat'),'table_','-v7.3');
elseif only_dur == 1
    save(strcat(exp_data_path,'GLM/training/',num2str(batdate),'_14650_datasetDistanceOnlyDur_allFlights_',num2str(nn),'.mat'),'table_','-v7.3');
else
    save(strcat(exp_data_path,'GLM/training/',num2str(batdate),'_14650_datasetComplete1_allFlights_',num2str(nn),'.mat'),'table_','-v7.3');
end

end
