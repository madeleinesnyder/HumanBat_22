function [globalCortexTime] = HumanBat_extractGlobalCortex(track_file,varargin)

%% Find global time for the cortex tracking data stream in b149f

% The 3s sync signal is on AnalogSignal channel 2 
% audiofile
% fsCortex = 120
%-----------------------------------------------------------------

%   *********USAGE EXAMPLES*****************

fsCortex = 120;
% Load tracking data
load(track_file);

% Find ttl audio peak indexes
[~,tracking_ttl_index] = findpeaks(AnalogSignals(:,2),'MinPeakHeight',0.3,'MinPeakDistance',fsCortex*2);
binaryTTLvector = NaN(1,length(AnalogSignals(:,2)));
binaryTTLvector(tracking_ttl_index)=0.36;
figure(); hold on; plot(AnalogSignals(:,2)); plot(binaryTTLvector,'*r');
trk_ttl_index_ds = tracking_ttl_index/fsCortex;
trk_ttl_index_us = trk_ttl_index_ds*1000000;
trk_ttl_time = [1:length(AnalogSignals(:,2))];
trk_ttl_time_ds = trk_ttl_time/fsCortex;
trk_ttl_us = trk_ttl_time_ds*1000000;

% Use local2GlobalTime to convert timestream of audio. 
[globalCortexTime] = local2GlobalTime(trk_ttl_index_us,trk_ttl_us);

end