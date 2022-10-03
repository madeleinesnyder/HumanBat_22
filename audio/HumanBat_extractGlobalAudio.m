function [globalMicTime] = HumanBat_extractGlobalAudio(audioConCat,ttlConCat,varargin)

%% Find global time for the audio data stream in b149f

% The sync signal is driven by a TTL on channel 6 of the recbuf in the  
% audiofile
% fsAudio = 192000 MOTU
% fsCortex = 120
%-----------------------------------------------------------------

%   *********USAGE EXAMPLES*****************

% Find ttl audio peak indexes
fsAudio=192000;
[~,audio_ttl_index] = findpeaks(ttlConCat,'MinPeakHeight',0.3,'MinPeakDistance',fsAudio*2);
binaryTTLvector = NaN(1,length(ttlConCat));
binaryTTLvector(audio_ttl_index)=0.36;
figure(); hold on; plot(ttlConCat); plot(binaryTTLvector,'*r');
audio_ttl_index_ds = audio_ttl_index/fsAudio;
audio_ttl_index_us = audio_ttl_index_ds*1000000;
audio_ttl_time = [1:length(ttlConCat)];
audio_ttl_time_ds = audio_ttl_time/fsAudio;
audio_ttl_us = audio_ttl_time_ds*1000000;

% Use local2GlobalTime to convert timestream of audio. 
[globalMicTime] = local2GlobalTime(audio_ttl_index_us,audio_ttl_us);
