function [B_audio_data] = HumanBat_alignEphystoCortexTTL(exp_data_path,audioConCat,ttlConCat)

    % Function to align the ephys data to the cortex TTL stream
    audio_things = dir(fullfile(strcat(exp_data_path,'/audio/audioConCat*')));
    ttl_things = dir(fullfile(strcat(exp_data_path,'/audio/ttlConCat*')));
    if length(audio_things) > 1
        disp("More than one audio trace extracted from this day")
    else
        load(strcat(audio_things.folder,'/',audio_things.name));
    end
    if length(ttl_things) > 1
        disp("More than one ttl trace extracted from this day")
    else
        load(strcat(ttl_things.folder,'/',ttl_things.name));
    end
    
    audio_fs = 192000;
    % Find the usec timestamp at which the first TTL in the audio stream occurred.
    [audio_ttl_peak_vals,audio_ttl_idxs] = findpeaks(ttlConCat(1:10000000),'MinPeakHeight',0.3,'MinPeakDistance',500000);
    audio_idx_vec = zeros(10000000,1); audio_idx_vec(audio_ttl_idxs) = 0.4;
    figure(); hold on; plot(ttlConCat(1:10000000)); plot(audio_idx_vec,'*r');
    first_audio_ttl = audio_ttl_idxs(1);
    first_audio_ttl_usec = first_audio_ttl/audio_fs*1e6;
    
    % Find the usec timestamp at which the last TTL in the audio stream occurred.
    clear audio_ttl_idxs audio_idx_vec 
    [audio_ttl_peak_vals,audio_ttl_idxs] = findpeaks(ttlConCat(end-10000000:end),'MinPeakHeight',0.3,'MinPeakDistance',500000);
    audio_idx_vec = zeros(10000000,1); audio_idx_vec(audio_ttl_idxs) = 0.4;
    figure(); hold on; plot(ttlConCat(end-10000000:end)); plot(audio_idx_vec,'*r');
    last_audio_ttl = size(ttlConCat,1) - 10000000 + audio_ttl_idxs(end);
    last_audio_ttl_usec = last_audio_ttl/audio_fs*1e6;
    
    clear ttlConCat
    %B_ttl_data = ttlConCat(first_audio_ttl:last_audio_ttl);
    B_audio_data = audioConCat(first_audio_ttl:last_audio_ttl);
    save(strcat(exp_data_path,'audio/B_audio_data_aligned.mat'),'B_audio_data','-v7.3');
end

