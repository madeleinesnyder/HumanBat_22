function [tot_dist,longest_flight,shortest_flight] = HumanBat_totalDistanceFlown(out)

    % Calculate total distance of bat flown during the session 
    tot_dist = 0;
    scale_factor = 1000;
    FD = [];
    for i=1:length(out.fstartxyz)
        % Make vector with all points in that flight
        flight_vector = out.trajectories_continuous(:,out.flight_starts_indx(i):out.flight_ends_indx(i));
        flight_vector_nonan = flight_vector(:,~all(isnan(flight_vector)));

        % Calculate distance between every point in that flight and sum
        flight_vector_dist = 0;
        for j=1:length(flight_vector_nonan)-1
            temp_dist = norm(flight_vector_nonan(:,j) - flight_vector_nonan(:,j+1))/scale_factor;
            flight_vector_dist = flight_vector_dist + temp_dist;
        end
        % Save as total distance
        FD(i) = flight_vector_dist;
        tot_dist = tot_dist + flight_vector_dist;
    end
    longest_flight = max(FD);
    shortest_flight = min(FD);
end