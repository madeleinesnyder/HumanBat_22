% Script to load in flights and trim the start and end times by hand

% Input:
    % cortex_flight_struct_resort{i}
    % cortex_r

% Output:
    % cortex_flight_struct_resort{i}.fstart_trimmed 
    % cortex_flight_struct_resort{i}.fend_trimmed

% MCS 6/22/22
% ============================================

% Specify batdate and cortex or cortex to get trimming
clear all;
batdate = 220418; exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
tracker = "cortex";
if tracker == "cortex"
    disp("Loading Pruned, Resorted Cortex struct and raw cortex data");
    load(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'));
    load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
    load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
end
close all;

%% Cortex case. Plot a flight and click on the start and end points
for i=1:length(cortex_flight_struct_resort)
    figure();  hold on; title(strcat("Spatial Profile of flight ",num2str(i))); xlim([-2900 2900]); ylim([-2300 2300]); zlim([0 2300])
    plot3(cortex_flight_struct_resort{i}.pos(1,:), cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:),'Color',[0.4 0.9 0]);
    scatter3(cortex_flight_struct_resort{i}.pos(1,1),cortex_flight_struct_resort{i}.pos(2,1),cortex_flight_struct_resort{i}.pos(3,1),'black','filled');
    
    p = cortex_flight_struct_resort{i}.pos;
    time = [1:length(cortex_flight_struct_resort{i}.pos)];
    v = zeros(length(time)-1,1);
    for j=1:length(time)-1
        v(j) = (p(j+1)-p(i))/(time(j+1)-time(j));
    end
    v=diff(p');
    figure(); hold on; title(strcat("Velocity profile of flight ",num2str(i)));
    plot(v);
    
    % Plot extensions of start and end
    start_lag = 200;
    stop_lag = 100;
    figure(); hold on; title(strcat("Extended Spatial profile of flight ",num2str(i)));
    plot3(cortex_flights.trajectoriesContinous(1,cortex_flight_struct_resort{i}.fstart-start_lag:cortex_flight_struct_resort{i}.fend+stop_lag),...
          cortex_flights.trajectoriesContinous(2,cortex_flight_struct_resort{i}.fstart-start_lag:cortex_flight_struct_resort{i}.fend+stop_lag),...
          cortex_flights.trajectoriesContinous(3,cortex_flight_struct_resort{i}.fstart-start_lag:cortex_flight_struct_resort{i}.fend+stop_lag),'Color',[0.3 0 0]);
    
    plot3(cortex_flights.trajectoriesContinous(1,cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend),...
          cortex_flights.trajectoriesContinous(2,cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend),...
          cortex_flights.trajectoriesContinous(3,cortex_flight_struct_resort{i}.fstart:cortex_flight_struct_resort{i}.fend),'Color',[0.8 0.8 0.8]);
    
    p_ = [cortex_flights.trajectoriesContinous(:,cortex_flight_struct_resort{i}.fstart-start_lag:cortex_flight_struct_resort{i}.fend+stop_lag)];
    time_ = [1:length(p_)];
    v_ = diff(p_');
    v_sm = smoothdata(v_,'gaussian',30); 
    f_ = figure(); hold on;
    title(strcat("Extended Flight ",num2str(i),". Click beginning and end, press Enter to trim."));
    plot(v_sm); 
    yline(0);
    hold off;
    
    % Get points from the figure and trim the flight to those points
    [start_trim,~] = getpts(f_);
    if ~isempty(start_trim)
        fstart_lagged = cortex_flight_struct_resort{i}.fstart-start_lag;    fend_lagged = cortex_flight_struct_resort{i}.fend+stop_lag;
        cortex_flight_struct_resort{i}.fstart_trimmed = floor(fstart_lagged+start_trim(1));
        cortex_flight_struct_resort{i}.fend_trimmed = ceil(fstart_lagged+start_trim(2));
        close all;
    else
        disp(strcat("No trimming flight ",num2str(i)));
        cortex_flight_struct_resort{i}.fstart_trimmed = cortex_flight_struct_resort{i}.fstart;
        cortex_flight_struct_resort{i}.fend_trimmed = cortex_flight_struct_resort{i}.fend;
    end

    p__= [cortex_flights.trajectoriesContinous(:,cortex_flight_struct_resort{i}.fstart_trimmed:cortex_flight_struct_resort{i}.fend_trimmed)];
    time__ = [1:length(p__)];
    v__ = diff(p__'); v__sm = smoothdata(v__,'gaussian',30); 
    f__ = figure(); hold on;
    title(strcat("Trimmed Flight ",num2str(i)));
    plot(v__sm);
    hold off;
    if mod(i,10) == 0
        save(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');
    end
end
save(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');

% OPTION 2 
for i=1:length(cortex_flight_struct_resort)
    cortex_flight_struct_resort{i}.fstart_trimmed = cortex_flight_struct_resort{i}.fstart;
    cortex_flight_struct_resort{i}.fend_trimmed = cortex_flight_struct_resort{i}.fend;
end
save(strcat(exp_data_path,'cortex/pruned_resorted_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');

