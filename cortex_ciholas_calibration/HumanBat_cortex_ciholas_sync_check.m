function [ciholas_cortex_synced] = HumanBat_cortex_ciholas_sync_check()

    % Function to check synchrony of ciholas and cortex on Blondie-only day
    % (This day should already be processed!)
    % Clear everything!
    clear all; 
    batdate = 220422;
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/14650/processed/',num2str(batdate),'/'); 
    % Load in the session that is Blondie only
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));
    load(strcat(exp_data_path,'cortex/aligned_bat_position_data.mat'));
    load(strcat(exp_data_path,'cortex/clustered_cortex_flights.mat'));
    load(strcat(exp_data_path,'ciholas/clustered_ciholas_flights.mat'));

    % Plot the data in 3space
    figure(); hold on;
    plot3(ciholas_r(:,1),ciholas_r(:,2),ciholas_r(:,3));
    plot3(cortex_flights.trajectoriesContinous(1,:)*1000,cortex_flights.trajectoriesContinous(2,:)*1000,cortex_flights.trajectoriesContinous(3,:)*1000);

    % Plot x coordinates
    figure(); hold on;
    plot(ciholas_r(:,1));
    plot(cortex_flights.trajectoriesContinous(1,:)*1000);

    % Plot the flights flight by flight
    n_cortex_flights = length(cortex_flights.id); n_ciholas_flights = size(ciholas_flights.flights,1);

    % From Cortex perspective
    for i=1:n_cortex_flights
        figure(); hold on;
        plot3(cortex_flights.trajectoriesContinous(1,cortex_flights.flight_starts_idx(i):cortex_flights.flight_ends_idx(i))*1000,...
        cortex_flights.trajectoriesContinous(2,cortex_flights.flight_starts_idx(i):cortex_flights.flight_ends_idx(i))*1000,...
        cortex_flights.trajectoriesContinous(3,cortex_flights.flight_starts_idx(i):cortex_flights.flight_ends_idx(i))*1000);
        plot3(ciholas_r(cortex_flights.flight_starts_idx(i):cortex_flights.flight_ends_idx(i),1),...
        ciholas_r(cortex_flights.flight_starts_idx(i):cortex_flights.flight_ends_idx(i),2),...
        ciholas_r(cortex_flights.flight_starts_idx(i):cortex_flights.flight_ends_idx(i),3)); hold off;    
    end

    % From ciholas perspective
    for i=1:n_ciholas_flights
        figure(); hold on; plot3(ciholas_r(ciholas_flights.flights{i,'smp1'}:ciholas_flights.flights{i,'smp2'},1),ciholas_r(ciholas_flights.flights{i,'smp1'}:ciholas_flights.flights{i,'smp2'},2),ciholas_r(ciholas_flights.flights{i,'smp1'}:ciholas_flights.flights{i,'smp2'},3));
        plot3(cortex_flights.trajectoriesContinous(1,ciholas_flights.flights{i,'smp1'}:ciholas_flights.flights{i,'smp2'})*1000,...
        cortex_flights.trajectoriesContinous(2,ciholas_flights.flights{i,'smp1'}:ciholas_flights.flights{i,'smp2'})*1000,...
        cortex_flights.trajectoriesContinous(3,ciholas_flights.flights{i,'smp1'}:ciholas_flights.flights{i,'smp2'})*1000);
        hold off;
    end

    prompt = "Does this plot look aligned? (Y) or (N)";
    g = input(prompt,'s');
    if g=="Y"
        ciholas_cortex_synced=1;
    else
        ciholas_cortex_synced=0;
    end
end

