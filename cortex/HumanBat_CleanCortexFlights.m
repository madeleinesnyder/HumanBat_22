% Function to clean up the Cortex Flights and get a continuous data stream
% equivilant to ciholas_r

load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
load(strcat(exp_data_path,'cortex/cortex_final_flight_structure.mat'));

% Re-oreder cortex flights
cv_f = [];
for i=1:length(cortex_flight_struct)
    cv_f = [cv_f,cortex_flight_struct{i}.fstart];
end
[im,imk] = sort(cv_f);
for i=1:length(cortex_flight_struct)
    cortex_flight_struct_resort{i} = cortex_flight_struct{imk(i)};
end

% Take a look at the differences bewteen cortex_flights.trajectoriesRaw and
% Continuous to see if the interpolation sucks
figure(); scatter3(cortex_flights.trajectoriesContinous(1,im(1):im(2)),cortex_flights.trajectoriesContinous(2,im(1):im(2)),cortex_flights.trajectoriesContinous(3,im(1):im(2)))
figure(); scatter3(cortex_flights.trajectoriesRaw(1,im(1):im(2)),cortex_flights.trajectoriesRaw(2,im(1):im(2)),cortex_flights.trajectoriesRaw(3,im(1):im(2)))

%% === Do the flights look ok? Prune out the shitty ones for now and split the ones that need splitting
to_split=[];
for i=1:length(cortex_flight_struct_resort)
    if isempty(cortex_flight_struct_resort{i})
        continue;
    else
        figure(); scatter3(cortex_flight_struct_resort{i}.pos(1,:),...
            cortex_flight_struct_resort{i}.pos(2,:),...
            cortex_flight_struct_resort{i}.pos(3,:));
        prompt = "Keep (K) or Toss (T) or Split (S)?";
        x = input(prompt,"s");
        if x=='T'
            cortex_flight_struct_resort{i} = [];
            close all;
        elseif x=='S'
            to_split = [to_split,i];
        else
            disp("Good flight, look at raster");
            close all;
        end
    end
end
save(strcat(exp_data_path,'cortex/pruned_cortex_bat_final_flight_structure.mat'),'cortex_flight_struct_resort','-v7.3');

%% split the flights that need splitting 
cortex_flight_struct_resort_split = cortex_flight_struct_resort;
start_lag = 100;
stop_lag = 60;
for i=1:length(to_split)
    p_ = [cortex_flights.trajectoriesContinous(:,cortex_flight_struct_resort{to_split(i)}.fstart-start_lag:cortex_flight_struct_resort{to_split(i)}.fend+stop_lag)];
    p_= p_';
    time_ = [1:length(p_)];
    v_ = diff(p_);
    f_ = figure(); hold on;
    title(strcat("Extended Flight ",num2str(to_split(i)),". Click beginning and end of FIRST flight, press Enter to trim."));
    plot(v_);
    yline(0);
    hold off;

    figure(); scatter3(cortex_flight_struct_resort{to_split(i)}.pos(1,:),...
            cortex_flight_struct_resort{to_split(i)}.pos(2,:),...
            cortex_flight_struct_resort{to_split(i)}.pos(3,:));
    
    % Get points from the figure and trim the flight to those points
    % (EXAMINE THIS CODE)
    [start_trim_1,~] = getpts(f_);
    fstart_lagged = cortex_flight_struct_resort{to_split(i)}.fstart-start_lag;    fend_lagged = cortex_flight_struct_resort{to_split(i)}.fend+stop_lag;
    cortex_flight_struct_resort_split{to_split(i)}.fstart_trimmed = floor(fstart_lagged+start_trim_1(1));
    cortex_flight_struct_resort_split{to_split(i)}.fend_trimmed = ceil(fstart_lagged+start_trim_1(2));
    close all;

    p_ = [cortex_flights.trajectoriesContinous(:,cortex_flight_struct_resort{to_split(i)}.fstart-start_lag:cortex_flight_struct_resort{to_split(i)}.fend+stop_lag)];
    p_= p_';
    time_ = [1:length(p_)];
    v_ = diff(p_);
    f_ = figure(); hold on;
    title(strcat("Extended Flight ",num2str(to_split(i)),". Click beginning and end of SECOND flight, press Enter to trim."));
    plot(v_);
    yline(0);
    hold off;
    [start_trim_2,~] = getpts(f_);
    fstart_lagged_1 = cortex_flight_struct_resort{to_split(i)}.fstart-start_lag;    fend_lagged_1 = cortex_flight_struct_resort{to_split(i)}.fend+stop_lag;
    cortex_flight_struct_resort_split{length(cortex_flight_struct_resort)+i}.fstart_trimmed = floor(fstart_lagged+start_trim_2(1));
    cortex_flight_struct_resort_split{length(cortex_flight_struct_resort)+i}.fend_trimmed = ceil(fstart_lagged+start_trim_2(2));
    close all;

    figure(); hold on; scatter3(cortex_flight_struct_resort{to_split(i)}.pos(1,:),...
            cortex_flight_struct_resort{to_split(i)}.pos(2,:),...
            cortex_flight_struct_resort{to_split(i)}.pos(3,:));
    scatter3(cortex_flight_struct_resort{to_split(i)}.pos(1,round(start_trim_1(1))),...
            cortex_flight_struct_resort{to_split(i)}.pos(2,round(start_trim_1(1))),...
            cortex_flight_struct_resort{to_split(i)}.pos(3,round(start_trim_1(1))),70,'r','filled');
    scatter3(cortex_flight_struct_resort{to_split(i)}.pos(1,round(start_trim_1(2))),...
            cortex_flight_struct_resort{to_split(i)}.pos(2,round(start_trim_1(2))),...
            cortex_flight_struct_resort{to_split(i)}.pos(3,round(start_trim_1(2))),70,'r','filled');

    scatter3(cortex_flight_struct_resort{to_split(i)}.pos(1,round(start_trim_2(1))),...
            cortex_flight_struct_resort{to_split(i)}.pos(2,round(start_trim_2(1))),...
            cortex_flight_struct_resort{to_split(i)}.pos(3,round(start_trim_2(1))),70,'g','filled');
    scatter3(cortex_flight_struct_resort{to_split(i)}.pos(1,round(start_trim_2(2))),...
            cortex_flight_struct_resort{to_split(i)}.pos(2,round(start_trim_2(2))),...
            cortex_flight_struct_resort{to_split(i)}.pos(3,round(start_trim_2(2))),70,'g','filled');
    close all;
end

%%
for i=2:length(cortex_flight_struct_resort)
    fstart = cortex_flight_struct_resort{i}.fstart;
    fend = cortex_flight_struct_resort{i}.fend;

    % Get ending time of the previous flight
    if i==1
        fstart_extended = fstart;
        fend_extended = cortex_flight_struct_resort{i+1}.fstart;
    else
        fstart_extended = cortex_flight_struct_resort{i-1}.fend;
        fend_extended = cortex_flight_struct_resort{i+1}.fstart;
    end

    figure(); hold on; xlim([-2.900 2.900]); ylim([-2.600 2.600]); zlim([0 2.300]); 
    scatter3(cortex_flights.trajectoriesContinous(1,fstart:fend_extended),cortex_flights.trajectoriesContinous(2,fstart:fend_extended),cortex_flights.trajectoriesContinous(3,fstart:fend_extended));
    scatter3(cortex_flights.trajectoriesContinous(1,fstart),cortex_flights.trajectoriesContinous(2,fstart),cortex_flights.trajectoriesContinous(3,fstart),60,'red','filled');

    figure(); hold on; xlim([-2.900 2.900]); ylim([-2.600 2.600]); zlim([0 2.300]);
    scatter3(cortex_flights.trajectoriesContinous(1,fstart_extended:fend),cortex_flights.trajectoriesContinous(2,fstart_extended:fend),cortex_flights.trajectoriesContinous(3,fstart_extended:fend));
    scatter3(cortex_flights.trajectoriesContinous(1,fstart_extended),cortex_flights.trajectoriesContinous(2,fstart_extended),cortex_flights.trajectoriesContinous(3,fstart_extended),60,'red','filled');

    figure(); hold on; xlim([-2.900 2.900]); ylim([-2.600 2.600]); zlim([0 2.300]);
    scatter3(cortex_flights.trajectoriesContinous(1,fstart:fend),cortex_flights.trajectoriesContinous(2,fstart:fend),cortex_flights.trajectoriesContinous(3,fstart:fend));
    scatter3(cortex_flights.trajectoriesContinous(1,fstart),cortex_flights.trajectoriesContinous(2,fstart),cortex_flights.trajectoriesContinous(3,fstart),60,'red','filled');

    

    % Interpolate and correct for the NaNs when the bat is at rest before
    % fstart.
    for j=1:(fstart-fstart_extended)
        if isnan(cortex_flights.trajectoriesContinous(1,fstart_extended+i-1))
            % TO DO fill in datapoint with rest location
            cortex_r(fstart_extended+i-1,1) = cortex_flights.trajectoriesContinous(1,fstart_extended);
            cortex_r(fstart_extended+i-1,2) = cortex_flights.trajectoriesContinous(2,fstart_extended);
            cortex_r(fstart_extended+i-1,3) = cortex_flights.trajectoriesContinous(3,fstart_extended);
        end
    end

end