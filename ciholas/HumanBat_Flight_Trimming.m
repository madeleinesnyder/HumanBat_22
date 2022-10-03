% Script to load in flights and trim the start and end times by hand

% Input:
    % ciholas_flight_struct_resort{i}
    % ciholas_r

% Output:
    % ciholas_flight_struct_resort{i}.fstart_trimmed 
    % ciholas_flight_struct_resort{i}.fend_trimmed

% MCS 6/22/22
% ============================================

% Specify batdate and ciholas or cortex to get trimming
clear all;
batdate = 220414; exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
tracker = "ciholas";
if tracker == "ciholas"
    disp("Loading Pruned, Resorted Ciholas struct and raw ciholas data");
    load(strcat(exp_data_path,'ciholas/pruned_resorted_ciholas_bat_final_flight_structure.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
else
    disp("Loading Cortex struct");
    load(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'));
end

%% Ciholas case. Plot a flight and click on the start and end points
for i=1:length(ciholas_flight_struct_resort)
    figure();  hold on; title(strcat("Spatial Profile of flight ",num2str(i))); xlim([-2900 2900]); ylim([-2300 2300]); zlim([0 2300])
    plot3(ciholas_flight_struct_resort{i}.pos(:,1), ciholas_flight_struct_resort{i}.pos(:,2),ciholas_flight_struct_resort{i}.pos(:,3),'Color',[0.4 0.9 0]);
    scatter3(ciholas_flight_struct_resort{i}.pos(1,1),ciholas_flight_struct_resort{i}.pos(1,2),ciholas_flight_struct_resort{i}.pos(1,3),'black','filled');
    
    p = ciholas_flight_struct_resort{i}.pos;
    time = [1:length(ciholas_flight_struct_resort{i}.pos)];
    v = zeros(length(time)-1,1);
    for j=1:length(time)-1
        v(j) = (p(j+1)-p(i))/(time(j+1)-time(j));
    end
    v=diff(p);
    figure(); hold on; title(strcat("Velocity profile of flight ",num2str(i)));
    plot(v);
    
    % Plot extensions of start and end
    start_lag = 100;
    stop_lag = 60;
    figure(); hold on; title(strcat("Extended Spatial profile of flight ",num2str(i)));
    plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart-start_lag:ciholas_flight_struct_resort{i}.fend+stop_lag,1),...
          ciholas_r(ciholas_flight_struct_resort{i}.fstart-start_lag:ciholas_flight_struct_resort{i}.fend+stop_lag,2),...
          ciholas_r(ciholas_flight_struct_resort{i}.fstart-start_lag:ciholas_flight_struct_resort{i}.fend+stop_lag,3),'Color',[0.3 0 0]);
    
    plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,1),...
          ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,2),...
          ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,3),'Color',[0.8 0.8 0.8]);
    
    p_ = [ciholas_r(ciholas_flight_struct_resort{i}.fstart-start_lag:ciholas_flight_struct_resort{i}.fend+stop_lag,:)];
    time_ = [1:length(p_)];
    v_ = diff(p_);
    f_ = figure(); hold on;
    title(strcat("Extended Flight ",num2str(i),". Click beginning and end, press Enter to trim."));
    plot(v_);
    yline(0);
    hold off;
    
    % Get points from the figure and trim the flight to those points
    [start_trim,~] = getpts(f_);
    fstart_lagged = ciholas_flight_struct_resort{i}.fstart-start_lag;    fend_lagged = ciholas_flight_struct_resort{i}.fend+stop_lag;
    ciholas_flight_struct_resort{i}.fstart_trimmed = floor(fstart_lagged+start_trim(1));
    ciholas_flight_struct_resort{i}.fend_trimmed = ceil(fstart_lagged+start_trim(2));
    close all;

    p__= [ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,:)];
    time__ = [1:length(p__)];
    v__ = diff(p__);
    f__ = figure(); hold on;
    title(strcat("Trimmed Flight ",num2str(i)));
    plot(v__);
    hold off;
    if mod(i,10) == 0
        save(strcat(exp_data_path,'ciholas/pruned_resorted_ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct_resort','-v7.3');
    end
end
save(strcat(exp_data_path,'ciholas/pruned_resorted_ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct_resort','-v7.3');

