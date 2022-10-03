% Script to plot landing and takeoff locations and save to a struct
clear all; close all;
cortex_ = 0;
ciholas_ = 1;
batdate=220411;  
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');

if ciholas_ 
    load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
    
    % Load in the trimmed or untrimmed data
    if exist(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'))
        load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
    else
        load(strcat(exp_data_path,'ciholas/pruned_resorted_ciholas_bat_final_flight_structure.mat'));
    end
elseif cortex_
    %load(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_final_flight_structure.mat'));
    load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
    
    % Load in the trimmed or untrimmed data
    if exist(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'))
        load(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
    else
        load(strcat(exp_data_path,'cortex/pruned_resorted_cortex_final_flight_structure.mat'));
    end
end

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
cutoff=600;

% Breifly plot tripod locations
figure(); scatter3(tripods(:,1),tripods(:,2),tripods(:,3));

%% For CIHOLAS!
if ciholas_

% Skeleton code for determining where a flight landed and took off
for i=1:length(ciholas_flight_struct_resort)
    flight_ending_coords(i,:) = mean(ciholas_r(ciholas_flight_struct_resort{i}.fend_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed+5,:));
    flight_starting_coords(i,:) = mean(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed-5:ciholas_flight_struct_resort{i}.fstart_trimmed,:));    
end
%%
for i=1:length(flight_ending_coords)
    dist_from_tripods = pdist2(tripods, flight_ending_coords(i,:), 'euclidean');
    tripod_landing=0;
    %if flight_ending_coords(i,3) > 1400
    %    flight_landing(i) = 8;
    %    continue
    %end
    for j=1:length(dist_from_tripods)
        if dist_from_tripods(j) < cutoff
            flight_landing(i) = j;
            tripod_landing = tripod_landing+1;
        end
    end
    if tripod_landing==0
        flight_landing(i)=10;
    elseif tripod_landing>1
        flight_landing(i)=11;
    end
end

for i=1:length(flight_starting_coords)
    dist_from_tripods = pdist2(tripods, flight_starting_coords(i,:), 'euclidean');
    tripod_takeoff=0;
    %if flight_starting_coords(i,3) > 1400
    %    flight_takeoff(i) = 8;
    %    continue
    %end
    for j=1:length(dist_from_tripods)
        if dist_from_tripods(j) < cutoff
            flight_takeoff(i) = j;
            tripod_takeoff = tripod_takeoff+1;
        end
    end
    if tripod_takeoff==0
        flight_takeoff(i)=10;
    elseif tripod_takeoff>1
        flight_takeoff(i)=11;
    end
end

for i=1:length(ciholas_flight_struct_resort)
    ciholas_flight_struct_resort{i}.tripod_landing = flight_landing(i);
    ciholas_flight_struct_resort{i}.tripod_takeoff = flight_takeoff(i);
end

%% Plot the clustered takeoff locations
unique_takeoff_locations = unique(flight_takeoff);
takeoff_colors = jet(length(unique(flight_takeoff))); 
figure(); hold on; %numbs = [1,2,3,4,5,6,8,9]; numbs_c = ['r','g','b','c','y','k','m','g','b','y','b','r'];
title("Takeoff Locations")
scatter3(flight_starting_coords(:,1),flight_starting_coords(:,2),flight_starting_coords(:,3),1,[0.9 0.9 0.9]);
for i=1:length(flight_starting_coords)
    for j=1:length(unique(flight_takeoff))
        if flight_takeoff(i) == unique_takeoff_locations(j)
            %scatter3(flight_starting_coords(i,1),flight_starting_coords(i,2),flight_starting_coords(i,3),numbs_c(j),'filled');
            scatter3(flight_starting_coords(i,1),flight_starting_coords(i,2),flight_starting_coords(i,3),20,takeoff_colors(j,:),'filled');
        end
    end
end
hold off;

unique_landing_locations = unique(flight_landing);
landing_colors = jet(length(unique_landing_locations)); 
figure(); hold on; 
title("Landing Locations")
scatter3(flight_ending_coords(:,1),flight_ending_coords(:,2),flight_ending_coords(:,3),5,[0.9 0.9 0.9]);
for i=1:length(flight_ending_coords)
    for j=1:length(unique_landing_locations)
        if j==length(unique_landing_locations)
            figure(); 
            plot3(cortex_flight_struct_resort{j}.pos(1,:),cortex_flight_struct_resort{j}.pos(2,:),cortex_flight_struct_resort{j}.pos(3,:));
        end
        if flight_landing(i) == unique_landing_locations(j)
            %scatter3(flight_ending_coords(i,1),flight_ending_coords(i,2),flight_ending_coords(i,3),numbs_c(j),'filled');
            scatter3(flight_ending_coords(i,1),flight_ending_coords(i,2),flight_ending_coords(i,3),20,landing_colors(j,:),'filled');
        end
    end
end
hold off;

% Save the tripod landing varables in the struct
if exist(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'))
    save(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct_resort','-v7.3');
else
    save(strcat(exp_data_path,'ciholas/pruned_resorted_ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct_resort','-v7.3');
end

%% FOR CORTEX
elseif cortex_

    % Skeleton code for determining where a flight landed and took off
for i=1:length(cortex_flight_struct_resort)
    flight_ending_coords(i,:) = mean(cortex_flight_struct_resort{i}.pos(:,end-3:end)');
    flight_starting_coords(i,:) = mean(cortex_flight_struct_resort{i}.pos(:,1:3)');
end

%%
for i=1:length(flight_ending_coords)
    dist_from_tripods = pdist2(tripods, flight_ending_coords(i,:), 'euclidean');
    tripod_landing=0;
    %if flight_ending_coords(i,3) > 1400
    %    flight_landing(i) = 8;
    %    continue
    %end
    for j=1:length(dist_from_tripods)
        if dist_from_tripods(j) < cutoff
            flight_landing(i) = j;
            tripod_landing = tripod_landing+1;
        end
    end
    if tripod_landing==0
        flight_landing(i)=11;
    elseif tripod_landing>1
        flight_landing(i)=12;
    end
end

for i=1:length(flight_starting_coords)
    dist_from_tripods = pdist2(tripods, flight_starting_coords(i,:), 'euclidean');
    tripod_takeoff=0;
    %if flight_starting_coords(i,3) > 1400
    %    flight_takeoff(i) = 8;
    %    continue
    %end
    for j=1:length(dist_from_tripods)
        if dist_from_tripods(j) < cutoff
            flight_takeoff(i) = j;
            tripod_takeoff = tripod_takeoff+1;
        end
    end
    if tripod_takeoff==0
        flight_takeoff(i)=11;
    elseif tripod_takeoff>1
        flight_takeoff(i)=12;
    end
end

for i=1:length(cortex_flight_struct_resort)
    cortex_flight_struct_resort{i}.tripod_landing = flight_landing(i);
    cortex_flight_struct_resort{i}.tripod_takeoff = flight_takeoff(i);
end

%% Plot the clustered takeoff locations
unique_takeoff_locations = unique(flight_takeoff);
takeoff_colors = jet(length(unique(flight_takeoff))); 
figure(); hold on; %numbs = [1,2,3,4,5,6,8,9]; numbs_c = ['r','g','b','c','y','k','m','g','b','y','b','r'];
title("Takeoff Locations")
scatter3(flight_starting_coords(:,1),flight_starting_coords(:,2),flight_starting_coords(:,3),1,[0.9 0.9 0.9]);
for i=1:length(flight_starting_coords)
    for j=1:length(unique_takeoff_locations)
        if flight_takeoff(i) == unique_takeoff_locations(j)
            %scatter3(flight_starting_coords(i,1),flight_starting_coords(i,2),flight_starting_coords(i,3),numbs_c(j),'filled');
            scatter3(flight_starting_coords(i,1),flight_starting_coords(i,2),flight_starting_coords(i,3),20,takeoff_colors(j,:),'filled');        end
    end
end
hold off;

% Plot all the weird flights
for i=1:length(flight_takeoff)
    if flight_takeoff(i) == 11
        figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
        scatter3(tripods(:,1),tripods(:,2),tripods(:,3),50,'b','filled');
        plot3(cortex_flight_struct_resort{i}.pos(1,:),cortex_flight_struct_resort{i}.pos(2,:),cortex_flight_struct_resort{i}.pos(3,:));
        scatter3(cortex_flight_struct_resort{i}.pos(1,1),cortex_flight_struct_resort{i}.pos(2,1),cortex_flight_struct_resort{i}.pos(3,1),50,'r','filled');
        hold off;
    end
end

unique_landing_locations = unique(flight_landing);
landing_colors = jet(length(unique_landing_locations)); 
figure(); hold on; 
title("Landing Locations")
scatter3(flight_ending_coords(:,1),flight_ending_coords(:,2),flight_ending_coords(:,3),5,[0.9 0.9 0.9]);
for i=1:length(flight_ending_coords)
    for j=1:length(unique_landing_locations)
        if flight_landing(i) == unique_landing_locations(j)
            %scatter3(flight_ending_coords(i,1),flight_ending_coords(i,2),flight_ending_coords(i,3),numbs_c(j),'filled');
            scatter3(flight_ending_coords(i,1),flight_ending_coords(i,2),flight_ending_coords(i,3),20,landing_colors(j,:),'filled');
        end
    end
end
hold off;

% Save the tripod landing varables in the struct
if exist(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_final_flight_structure.mat'))
    save(strcat(exp_data_path,'cortex/pruned_resorted_trimmed_cortex_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');
else
    save(strcat(exp_data_path,'cortex/pruned_resorted_cortex_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');
end
end