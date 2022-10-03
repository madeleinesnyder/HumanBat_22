function [] = HumanBat_brainMovies(movie_path,output_filename,batdate,logger,unit,camera_angle, starting_time,ending_time)

% Called by HumanBat_shell_brainMovies

% Inputs:
    % movie_path: path to the short movie clip specified by the starting and ending times
    % output_filename: the output movie file
    % batdate (i.e. 220407)
    % logger (i.e. 13)
    % unit (i.e. 5)
    % camera angle ("wall","one","two","top","front")
    % start time in seconds (i.e. 490 (aka 8:10))
    % duration in seconds (i.e. 300 (aka 5 minutes))
% Outputs:
    % Writes and saves the avi that has audio and video

% MCS 6/22/22
% ============================================

exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/',num2str(batdate),'/');
% Load in the ephys data
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
% Load in the flight data
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data.mat'));

% Length of session in seconds
session_length = length(ciholas_r)/120;
session_length_usec = session_length*1e6;

% Get spike sounds binary vector
[sound_clip] = HumanBat_spikes2sound(B_ephys_data.TT_unit(unit).AlignedTimestamps,session_length_usec,starting_time,ending_time);

% Load in video 
% To split a new video: ffmpeg -i input.mp4  -ss starting time in minute mark (i.e. 00:08:10) -t amount of time in new clip (i.e. 00:05:00) -c:v copy -c:a copy output.mp4
videoFReader = vision.VideoFileReader(movie_path);
vv = VideoReader(movie_path);

%% Write the audio to the video
fps = vv.NumFrames/vv.Duration;
audiops = 100000;
samplesPerFrame = audiops/fps;
videoFWriter = vision.VideoFileWriter('Filename',output_filename,'AudioInputPort',true,'FrameRate',fps,'FileFormat','AVI');
framenum=0;
while ~isDone(videoFReader)
    disp(strcat(num2str(framenum),'/',num2str(vv.NumFrames)));
    Fframe = step(videoFReader);
    thus_audio = sound_clip.binary(framenum*samplesPerFrame + 1 : min(end, (framenum+1)*samplesPerFrame))';
    step(videoFWriter,Fframe,thus_audio);
    framenum=framenum+1;
end
release(videoFReader);
release(videoFWriter);
end


