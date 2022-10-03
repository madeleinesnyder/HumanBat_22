% Extract c3d files (either pipeline or in directory) with c3d_batchextract.m


% Format the tracking Marker data
[Location, Location_interp] = ImBat_formatTracking(Markers);
Location_times = [1:length(AnalogSignals)];

% Segment the flights
[out] =  ImBat_SegTrajectories(Location,Location_times);

% Plot the flights to see by eye and check
num_flights = length(out.flight_starts_indx); figure();
disp(strcat("14592 performed "," ",num2str(num_flights)," ","during flightroom b149f task."));
for i=1:num_flights
    subplot(round(sqrt(num_flights))+1,round(sqrt(num_flights))+1,i); hold on;
    scatter3(out.trajectories_continuous(1,out.flight_starts_indx(i):out.flight_ends_indx(i)),out.trajectories_continuous(2,out.flight_starts_indx(i):out.flight_ends_indx(i)),out.trajectories_continuous(3,out.flight_starts_indx(i):out.flight_ends_indx(i)));
end

% Plot all trajectories flown during the task
figure(); scatter3(out.trajectories_continuous(1,:),out.trajectories_continuous(2,:),out.trajectories_continuous(3,:))

% Plot total distance flown during the task
[total_distance_flown,longest_flight,shortest_flight] = HumanBat_totalDistanceFlown(out);
disp(strcat("Total Distance Flown:"," ",num2str(total_distance_flown)," ","meters."));
disp(strcat("Longest Flight:"," ",num2str(longest_flight)," ","meters."));
disp(strcat("Shortest Flight:"," ",num2str(shortest_flight)," ","meters."))

% Plot a timeline of when the flights were flown
binaryFlightVector = zeros(1,length(AnalogSignals));
binaryFlightVector(out.flight_starts_indx)=1;
figure(); stem(binaryFlightVector); title("Timeline of Flights during session");

% Sort the flights
flightPaths = HumanBat_ClusterFlights(out,AnalogSignals);
