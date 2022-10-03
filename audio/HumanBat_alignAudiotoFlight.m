function HumanBat_alignAudiotoFlight(DATE)

% Load in flight data 
track_file = strcat(pwd,'/data/processed/',DATE,'/b149f/cortex/',DATE,'_14592_tracking_1_track.mat');
load(track_file);

% Get global cortex time if needed 
[globalCortexTime] = HumanBat_extractGlobalCortex(track_file);

% Load in audio data (ttl and mic)
audio_dir = dir(strcat(pwd,'/data/raw/',DATE,'/b149f/audio'));
mic_file = strcat(pwd,'/data/raw/',DATE,'/b149f/audio/',audio_dir(3).name,'/audioConCat_1.mat');
audio_ttl_file = strcat(pwd,'/data/raw/',DATE,'/b149f/audio/',audio_dir(3).name,'/ttlConCat.mat');
load(mic_file);
load(audio_ttl_file);

% Get global mic and audio_ttl time if needed
[globalMicTime] = HumanBat_extractGlobalAudio(audioConCat,ttlConCat);

% Plot aligning things

% TODO
seconds_to_plot = 300;
audio_samples = 100*192000;
cortex_samples = 100*120;

% Plot analog Signals and align to first TTL
disp(strcat('Cortex time is:'," ",num2str(length(AnalogSignals)/fsCortex/60)," ",'minutes'));
disp(strcat('Audio time is:'," ",num2str(length(audioConCat)/fsAudio/60)," ",'minutes'));
figure(); 
subplot(3,1,1); hold on; plot(AnalogSignals(1:cortex_samples,2)); title(strcat('Cortex ttl first '," ",num2str(seconds_to_plot)," ","seconds")); hold off;
subplot(3,1,2); hold on; plot(ttlConCat(1:audio_samples)); title(strcat('TTL Audio first '," ",num2str(seconds_to_plot)," ","seconds")); hold off;
subplot(3,1,3); hold on; plot(audioConCat(1:audio_samples)); title(strcat('Mic Audio first '," ",num2str(seconds_to_plot)," ","seconds")); hold off;



end