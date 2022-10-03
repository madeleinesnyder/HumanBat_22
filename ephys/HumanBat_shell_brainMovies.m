function [] = HumanBat_shell_brainMovies(batdate,logger,unit,camera_angle,start_time,dur,sanity_check_flag)

% Creates movie with spikes overlaid like Hubel&Weisel cat bar videos
    % Executes terminal commands to split the movies as you wish
    % Extracts spikes and binarizes them as audio at 100000Hz
    % Writes file as avi
    % Compresses to mp4

% Inputs:
    % batdate (i.e. 220407)
    % logger (i.e. 13)
    % unit (i.e. 5)
    % camera angle ("wall","one","two","top","front")
    % start time in seconds (i.e. 490 (aka 8:10))
    % duration in seconds (i.e. 300 (aka 5 minutes))
    % sanity_check_flag (if 1, plot the times of the spikes, if 0, dont')
% Outputs:
    % .mp4 file of aligned video and 

% Sample run:
    % HumanBat_shell_brainMovies(220407,13,5,"top",490,180);

% MCS 6/22/22
% ============================================

% Dummy inputs
% batdate=220407; logger=13; unit=5; camera_angle="front"; start_time=960; dur=180;

% Load paths
addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/HumanBat/ephys'));
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
exp_movie_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/raw/',num2str(batdate),'/cameras/',camera_angle);
video_dir = dir(fullfile(strcat(exp_movie_path,'/Basler*')));
video_file = video_dir.name;
output_video_file = strcat(video_dir.folder,'/unit_',num2str(unit),'_t_',num2str(start_time),'_',num2str(start_time+dur),'.mp4');
output_video_file_audiocombined = strcat(video_dir.folder,'/unit_',num2str(unit),'_t_',num2str(start_time),'_',num2str(start_time+dur),'_withAudio.avi');
output_video_file_compressed = strcat(video_dir.folder,'/unit_',num2str(unit),'_t_',num2str(start_time),'_',num2str(start_time+dur),'_withAudio.mp4');

% ========================================================

% 1. Split desired video
start_time_formatted = duration(0,0,start_time);
dur_formatted = duration(0,0,dur);
command_split = strcat("ffmpeg -i ", strcat(video_dir.folder,'/',video_file), " -ss ",string(start_time_formatted)," -t ",string(dur_formatted)," -c:v copy -c:a copy ", output_video_file);
outpt_split = system(command_split);

% 2. Run script to combine audio and video
HumanBat_brainMovies(output_video_file,output_video_file_audiocombined,batdate,logger,unit,camera_angle,start_time,start_time+dur);

% 3. Compress video clip
command_compress = strcat("ffmpeg -i ",output_video_file_audiocombined," -c:v libx265 -x265-params lossless=1 -c:a libfdk_aac -b:a 128k -y ",output_video_file_compressed);
outpt_compress = system(command_compress);

% 4. Delete .AVI file
avi_files = dir(fullfile(strcat(exp_movie_path,'/*.avi')));
avi_files_to_delete = strcat(avi_files.folder,'/',avi_files.name);
delete(avi_files_to_delete)

% 5. Sanity check that the raster lines up with the auditory spikes
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
start_time_usec = start_time*1e6; end_time_usec = (start_time+dur)*1e6;
ephys_time_vec = zeros(1,end_time_usec-start_time_usec);
ephys_timestamps_in_time_vec = B_ephys_data.TT_unit(unit).AlignedTimestamps(B_ephys_data.TT_unit(unit).AlignedTimestamps > start_time_usec);
ephys_timestamps_in_time_vec = ephys_timestamps_in_time_vec(ephys_timestamps_in_time_vec < end_time_usec);
ephys_timestamps_in_time_vec = ephys_timestamps_in_time_vec - start_time_usec;
ephys_time_vec(round(ephys_timestamps_in_time_vec)) = 1;
% Plot 
figure(); scatter([1:length(ephys_time_vec)],ephys_time_vec); title(strcat("Unit ",num2str(unit)));

end