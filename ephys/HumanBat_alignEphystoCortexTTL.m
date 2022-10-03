function [B_ephys_data] = HumanBat_alignEphystoCortexTTL(B_ephys_data,cortex_ttl_seconds,ephys_TTL)

% Function to align the ephys data to the cortex TTL stream

% Find the TTL at which the first ttl was detectd in the ephys stream
ephys_ttl_idxs = find(contains(ephys_TTL.event_types_and_details, 'Digital in. Digital input port status = '));%0x29ef0401; Digital input event status = 0x00000001; Digital in rising edge on pin number 1. Event name: Pin1Event;'));
ephys_TTL_seconds = ephys_TTL.event_timestamps_usec(ephys_ttl_idxs(1));
ephys_last_ttl_idx = ephys_TTL.event_timestamps_usec(ephys_ttl_idxs(end));

% All timestamps in the ephys data stream are in 1e6 seconds. They are the
% 1e6 second at which that cell (i) fired.
% First scale 1e6 seconds to 0 in the ephys
for i=1:length(B_ephys_data.TT_unit)
    B_ephys_data.TT_unit(i).AlignedTimestamps = B_ephys_data.TT_unit(i).Timestamps-ephys_TTL_seconds;
    B_ephys_data.TT_unit(i).WRONG = B_ephys_data.TT_unit(i).Timestamps-B_ephys_data.TT_unit(i).Timestamps(1);
end

end

