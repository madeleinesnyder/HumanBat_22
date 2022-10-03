function [] = HumanBat_SI_red_dots(batdate,logger,unit,ciholas_flight_struct_resort,rdc)
% Script to plot the spatial infomration in a flight-centric way (Red
% dots plots)

% Inputs:
    % unit: putative unit to examine SI of (i.e. 5)
    % batdate: date of session
    % logger: which logger data to examine
% Outputs:
    % Plot of all flights with spikes overlaid in red
    % Plot of all clusterable flights with spikes overlaid in red
    % Plots of the top 5 spike-y flights 

% MCS 6/22/22
% ============================================

%close all;
%unit=13; 
%batdate=220408;  logger=13;
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));

% Load in the ciholas flight structure that has been pruned of the shitty
% flights and resorted according to cluster identity
%load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure.mat'));
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));

%% For every flight, plot the red dots on that flight
figure(); hold on; title("Flight with ephys plotted on top");
for i=1:length(ciholas_flight_struct_resort)
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
hold off;

%% For every flight in cluster clus plot with red dots

figure(); hold on; title(strcat("Flight clusters 2:end with ephys plotted on top")); xlim([-2900 2900]); ylim([-2500 2500]); zlim([0 2300]);
for p = 2:ciholas_flight_struct_resort{end}.fclus
    clus=p;
    for i=1:length(ciholas_flight_struct_resort)
        if ciholas_flight_struct_resort{i}.fclus==clus
            plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,1),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,2),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,3),'Color','b');
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
        else
            continue
        end
    end
end

%% For every flight in cluster clus plot with red dots INDIVIDUAL FIGURES

for p = 1:length(rdc)
    clus = rdc(p);
    figure(); hold on; title(strcat("Flight cluster ",num2str(clus))); xlim([-2900 2900]); ylim([-2500 2500]); zlim([0 2300]);
    for i=1:length(ciholas_flight_struct_resort)
        if ciholas_flight_struct_resort{i}.fclus==clus
            plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,1),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,2),...
                  ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed:ciholas_flight_struct_resort{i}.fend_trimmed,3),'Color',[0.8 0.8 0]);
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
        else
            continue
        end
    end
    hold off;
end


%% Find the flights with the most firing
for i=1:length(ciholas_flight_struct_resort)
    num_spikes_per_flight(i) = length(ciholas_flight_struct_resort{i}.ephys_trimmed{unit});
    num_spikes_per_flight_norm_ts(i) = length(ciholas_flight_struct_resort{i}.ephys_trimmed{unit})/(length(ciholas_flight_struct_resort{i}.pos)/120);
end

[spikiest,spikiest_idx] = sort(num_spikes_per_flight);  [spikiest_norm,spikiest_idx_norm] = sort(num_spikes_per_flight_norm_ts);

% Plot the spikiest flights
spikiest_flights = spikiest_idx_norm(end-5:end);
for j=1:length(spikiest_flights)
    figure(); hold on; title("Spikiest flights; K yellow M blue"); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    i = spikiest_flights(j);
    %plot3(ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,1),...
    %                   ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,2),...
    %                   ciholas_r(ciholas_flight_struct_resort{i}.fstart:ciholas_flight_struct_resort{i}.fend,3),'Color',[0.8 0 0]);
    plot3(ciholas_flight_struct_resort{i}.pos(:,1), ciholas_flight_struct_resort{i}.pos(:,2),ciholas_flight_struct_resort{i}.pos(:,3),'Color',[0.8 0.9 0]);
    % Transform ephys timestamps to ciholas time
    % Make artificial vector the time duration of ciholas
    flight_dur = (ciholas_flight_struct_resort{i}.fend_trimmed -  ciholas_flight_struct_resort{i}.fstart_trimmed);
    flight_vec = [1:flight_dur]./120;
    if ciholas_flight_struct_resort{i}.ephys_trimmed{unit} == 0
        disp("No unit firing on this flight");
        continue;
    else
        for jj=1:length(ciholas_flight_struct_resort{i}.ephys_trimmed{unit})
            nearest_ephys_point = dsearchn(flight_vec',ciholas_flight_struct_resort{i}.ephys_trimmed{unit}(jj)/1e6);
            try
                scatter3(ciholas_flight_struct_resort{i}.pos(nearest_ephys_point,1),ciholas_flight_struct_resort{i}.pos(nearest_ephys_point,2),ciholas_flight_struct_resort{i}.pos(nearest_ephys_point,3),'r','filled');
            catch
                disp("Trimmed flight is LONGER than original flight. Missing some points plotted");
            end
        end
    end
    % Plot human positions
    scatter3(ciholas_flight_struct_resort{i}.pos_K(:,1), ciholas_flight_struct_resort{i}.pos_K(:,2),ciholas_flight_struct_resort{i}.pos_K(:,3),'y');
    scatter3(ciholas_flight_struct_resort{i}.pos_M(:,1), ciholas_flight_struct_resort{i}.pos_M(:,2),ciholas_flight_struct_resort{i}.pos_M(:,3),'b');
    scatter3(ciholas_flight_struct_resort{i}.pos_cortex_bat(1,:),ciholas_flight_struct_resort{i}.pos_cortex_bat(2,:),ciholas_flight_struct_resort{i}.pos_cortex_bat(3,:),'green');
    scatter3(ciholas_flight_struct_resort{i}.pos(1,1),ciholas_flight_struct_resort{i}.pos(1,2),ciholas_flight_struct_resort{i}.pos(1,3),'black','filled');
end
hold off;

end
