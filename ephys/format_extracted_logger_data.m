function format_extracted_logger_data(logger_dir)
%format_extracted_logger_data Loads extracted data from individual .mat files and saves into a
%(num_channels, num_samples) format in both .mat and .bin
%Also adds global sample timestamps in usec referenced to first TTL (t = 0).

disp(logger_dir)
data = load_extracted_data_logger(fullfile(logger_dir, 'extracted_data'));


if(~isfile(fullfile(logger_dir, 'logger_data.mat')))
    disp("Saving logger data to .mat")

    % Putative local sample timestamps
    t0 = data.timestamps_first_samples_logger_usec{1,1}(1);
    local_sample_ts = linspace(t0, t0+(length(data.csc)-1)*data.sampling_period_usec, length(data.csc));
    % Convert local sample timestamps to global timestamps (relative to
    % first master9 ttl)
    data.global_sample_timestamps_usec = local2GlobalTime(data.ttl_timestamps_usec,local_sample_ts,'global_ttl_interval_us', 3e6);
    save(fullfile(logger_dir, 'logger_data.mat'), 'data','-v7.3')
end

if(~isfile(fullfile(logger_dir, 'logger_data.bin')))
    disp("Saving logger data to .bin")
    fileID = fopen(fullfile(logger_dir, 'logger_data.bin'), 'w');
    fwrite(fileID, data.csc.', 'int16');
    fclose(fileID)
end

%if(~isfile(fullfile(logger_dir, 'logger_data.mat')))
%    disp("Saving CAR and HP filtered logger data to .bin")
%    for i=1:16
%        disp(sprintf('highpass filtered channel %d',i));
%        data.csc(:,i) = highpass(double(data.csc(:, i)), 300, 31250);
%    end
%    med_csc = median(data.csc.',1);
%    fileID = fopen(fullfile(logger_dir, 'logger_data_car.bin'), 'w');
%    fwrite(fileID, data.csc.' - med_csc, 'int16');
%    fclose(fileID)
%end


end