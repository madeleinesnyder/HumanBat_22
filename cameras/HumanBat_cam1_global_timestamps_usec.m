function [global_cam1_usec] = HumanBat_cam1_global_timestamps_usec(path_to_mp4)

V = VideoReader(path_to_mp4);
D = V.duration; % Get length of recording in seconds
frames = V.NumFrames;

% Load in the ciholas data


end