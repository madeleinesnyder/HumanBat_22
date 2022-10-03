function [sound_clip] = HumanBat_spikes2sound(spike_timestamps,trial_length,starting_time,ending_time)

% Function to convert spikes to sound (binary vector)

% Input:
    % spike_timestamps: the millisecond time of a given unit's spikes, aligned to all other data streams (or, video data stream)
    % trial_length: length of the trial in milliseconds (or segment of the trial) from which the spike_timestamps have been extracted and aligned to 
    % starting_time: the start time of the clip 
    % ending_time: the end time of the clip 

% Output:
    % sound_clip struct with 
    % sound_clip.binary = [0 0 . . . 1 0]
    % sound_clip.rate = 100000

% MCS 6/22/22
% ============================================

max_matlab_rate = 100000;
real_sample_rate = 1e6;

% Truncate the spike timestamps vector to the desired length
trial_length_trunc = ending_time*real_sample_rate-starting_time*real_sample_rate;
spike_times_trunc = spike_timestamps(spike_timestamps > starting_time*real_sample_rate);
spike_times_trunc = spike_times_trunc(spike_times_trunc <  ending_time*real_sample_rate);

% Bring it to seconds then to max rate
trial_length_max = round(trial_length_trunc/real_sample_rate*max_matlab_rate);
for i=1:length(spike_times_trunc)
    spike_timestamps_max(i) = round(spike_times_trunc(i)/real_sample_rate*max_matlab_rate);
end

% Create binary vector of the rescaled trial length where 1's are spikes.
% This will be zerod for the truncated piece
spike_vec = zeros(1,trial_length_max);
spike_timestamps_max_zerod = spike_timestamps_max-starting_time*max_matlab_rate;
%spike_vec(1) = 1;
spike_vec(spike_timestamps_max_zerod) = 1;

% Play sound from time time 
sound_clip.binary = spike_vec;
sound_clip.rate = max_matlab_rate;
%sound(spike_vec,max_matlab_rate)

end