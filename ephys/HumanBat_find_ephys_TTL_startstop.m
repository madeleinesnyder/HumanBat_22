function [TTL_start_stop] = HumanBat_find_ephys_TTL_startstop(eventfile)

    % Function to make a matfile with start and stop of ttl for flightroom
    % session. This file will be used in NlxtMattSpike conversion of a
    % sorted spikefile to a struct with the spike features. 

    load(eventfile);
    TTL_timestamp_idxs = [];
    for i=1:length(event_types_and_details)
        temp_match = regexp(event_types_and_details{i},'Digital in rising edge on pin number 1');
        if isempty(temp_match)
            continue
        else
            TTL_timestamp_idxs = [TTL_timestamp_idxs,event_timestamps_usec(i)];
        end
    end

    % Check that the eventfile is valid and has all TTLs and they are all
    % evenly spaced
    figure(); subplot(2,1,1); plot(TTL_timestamp_idxs); subplot(2,1,2); plot(diff(TTL_timestamp_idxs));
    % Back of the envalope for amount of time
    disp(strcat("Session length according to phys TTLs:"," ",num2str((length(TTL_timestamp_idxs)*3)/60)," ","minutes"));
   
    % Return vector that contains first and last TTLs
    TTL_timestamps = [TTL_timestamp_idxs(1),TTL_timestamp_idxs(end)];
    save('TTL_timestamps','TTL_timestamps');
end
