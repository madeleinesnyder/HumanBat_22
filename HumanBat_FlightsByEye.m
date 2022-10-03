% Script to open and play the video from the correct time point for
% Mx_flights and Kx_flights

Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)

% Camera boundings 
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

% Cams are (1 - one, 2 - two, 3 - wall, 4 - corner, 5 - front, 6 - top)
% Cams map to tripods: (1 - 3, 2 - 5, 3 - 6, 4 - 1+2, 5 - 4+7, 6 - all)
cam1 = tripods(3,:); cam2 = tripods(5,:); cam3 = tripods(6,:); cam4 = mean([tripods(1,:);tripods(2,:)]); cam5 = tripods(4,:);
cams = [cam1;cam2;cam3;cam4;cam5];

cam_dict = {'one','two','wall','corner','front','top'};

%% CIHOLAS
clus = clus;
figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
for i=1:length(Flight_Group_Matrix.Clusters{clus})
    flightnum = Flight_Group_Matrix.Clusters{clus}(i);
    plot3(ciholas_flight_struct_resort{flightnum}.pos(:,1),ciholas_flight_struct_resort{flightnum}.pos(:,2),ciholas_flight_struct_resort{flightnum}.pos(:,3),'k');
    scatter3(ciholas_flight_struct_resort{flightnum}.pos(:,1),ciholas_flight_struct_resort{flightnum}.pos(:,2),ciholas_flight_struct_resort{flightnum}.pos(:,3),'k');
end
title(strcat("Cluster ",num2str(clus)));
hold off;

% Find out which cameras we need to look at for the takeoff and the landing
takeoff_camera_type = [];
landing_camera_type = [];
for i=1:length(Flight_Group_Matrix.Clusters{clus})
    flightnum = Flight_Group_Matrix.Clusters{clus}(i);
    landing_coords = [mean(ciholas_flight_struct_resort{flightnum}.pos(end-10:end,1)),mean(ciholas_flight_struct_resort{flightnum}.pos(end-10:end,2)),mean(ciholas_flight_struct_resort{flightnum}.pos(end-10:end,3))];
    cam_dist = sum(abs(cams-landing_coords),2);
    best_landing_cam = find(cam_dist == min(cam_dist));
    takeoff_coords = [mean(ciholas_flight_struct_resort{flightnum}.pos(1:10,1)),mean(ciholas_flight_struct_resort{flightnum}.pos(1:10,2)),mean(ciholas_flight_struct_resort{flightnum}.pos(1:10,3))];
    cam_dist = sum(abs(cams-takeoff_coords),2);
    best_takeoff_cam = find(cam_dist == min(cam_dist));
end

to_list_m = []; to_list_k = [];
for j=1:length(Mx_flights)
    idx = Mx_flights(j);
    to_list_m = [to_list_m,ciholas_flight_struct_resort{idx}.fstart_trimmed/120];
end

for j=1:length(Kx_flights)
    idx = Kx_flights(j);
    to_list_k = [to_list_k,ciholas_flight_struct_resort{idx}.fstart_trimmed/120];
end

% Add in coding for Kevin (2) versus Madelien (1) and sort the full list
to_list_m = to_list_m'/60;  to_list_m = [to_list_m,repmat(1,length(to_list_m),1),Mx_flights'];  
to_list_k = to_list_k'/60;  to_list_k = [to_list_k,repmat(2,length(to_list_k),1),Kx_flights'];   
to_list = [to_list_m;to_list_k];
[sorted_to_list,idx] = sortrows(to_list,1);

disp(strcat("Video: ",num2str(batdate)));
disp(strcat("Takeoff Camera View: ",cam_dict{best_takeoff_cam}));
disp(strcat("Landing Camera View: ",cam_dict{best_landing_cam}));
disp("Times: ");
for i=1:size(sorted_to_list,1)
    if sorted_to_list(i,2) == 1
        disp(strcat(num2str(sorted_to_list(i,1))," M"));
    else
        disp(strcat(num2str(sorted_to_list(i,1))," K"));
    end  
end




%% FOR CORTEX

% Find out which cameras we need to look at for the takeoff and the landing
takeoff_camera_type = [];
landing_camera_type = [];
for i=1:length(Flight_Group_Matrix.Clusters{clus})
    flightnum = Flight_Group_Matrix.Clusters{clus}(i);
    landing_coords = [mean(cortex_flight_struct_resort{flightnum}.pos(1,end-10:end)),mean(cortex_flight_struct_resort{flightnum}.pos(2,end-10:end)),mean(cortex_flight_struct_resort{flightnum}.pos(3,end-10:end))];
    cam_dist = sum(abs(cams-landing_coords),2);
    best_landing_cam = find(cam_dist == min(cam_dist));
    takeoff_coords = [mean(cortex_flight_struct_resort{flightnum}.pos(1,1:10)),mean(cortex_flight_struct_resort{flightnum}.pos(2,1:10)),mean(cortex_flight_struct_resort{flightnum}.pos(3,1:10))];
    cam_dist = sum(abs(cams-takeoff_coords),2);
    best_takeoff_cam = find(cam_dist == min(cam_dist));
end

to_list_m = []; to_list_k = [];
for j=1:length(Mx_flights)
    idx = Mx_flights(j);
    to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart/120];
end

for j=1:length(Kx_flights)
    idx = Kx_flights(j);
    to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart/120];
end

% Add in coding for Kevin (2) versus Madelien (1) and sort the full list
to_list_m = to_list_m'/60;  to_list_m = [to_list_m,repmat(1,length(to_list_m),1),Mx_flights'];  
to_list_k = to_list_k'/60;  to_list_k = [to_list_k,repmat(2,length(to_list_k),1),Kx_flights'];   
to_list = [to_list_m;to_list_k];
[sorted_to_list,idx] = sortrows(to_list,1);

disp(strcat("Video: ",num2str(batdate)));
disp(strcat("Takeoff Camera View: ",cam_dict{best_takeoff_cam}));
disp(strcat("Landing Camera View: ",cam_dict{best_landing_cam}));
disp("Times: ");
for i=1:size(sorted_to_list,1)
    if sorted_to_list(i,2) == 1
        disp(strcat(num2str(sorted_to_list(i,1))," M"));
    else
        disp(strcat(num2str(sorted_to_list(i,1))," K"));
    end  
end
