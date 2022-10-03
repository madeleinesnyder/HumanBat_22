extracted_cdp = load('ciholas/extracted_220104_cdp_1.mat');
x = extracted_cdp.tag_data_filt{1}(1:end,3);
y = extracted_cdp.tag_data_filt{1}(1:end,4);
cdp_fs = 100;


extracted_cortex = load('cortex/Generated_C3D_files/processed/220104_14705_tracking_1_track.mat');
cortex_pos = extracted_cortex.avgMarkerPos; % cortex tracked pos
cortex_fs = extracted_cortex.AnalogFrameRate; % cortex sampling rate
ciholas_fs = 100;

ciholas_first_sample_ind = find(extracted_cdp.sync_data(:,8) == 0);

ttl_3s = extracted_cortex.AnalogSignals(1:end, 2); % 3 sec ttl (ephys)
ttl_fib = extracted_cortex.AnalogSignals(1:end, 3); % fib ttl (ciholas)
[R,LT,UT,LL,UL] = risetime(ttl_3s,cortex_fs);
first_3s_ttl_ms = LT(1)*1000; % Timestamp of first 3s TTL in ms
last_3s_ttl_ms = LT(end)*1000; % Timestamp of last 3s TTL in ms
[R,LT,UT,LL,UL] = risetime(ttl_fib,cortex_fs);
first_fib_ttl_ms = LT(1)*1000; % Timestamp of first fib TTL in ms
last_fib_ttl_ms = LT(end)*1000; % Timestamp of last fib TTL in ms
first_ttl_sample_ind_120hz = round(first_fib_ttl_ms*cortex_fs/1e3) + round(1.05*cortex_fs); % Sample index at 120hz corresponding to first ttl 
%last_ttl_sample_ind_120hz = round(last_cortex_ttl_ms*cortex_fs/1e3); % Sample index at 120hz corresponding to first ttl 
cortex_ttl_ms = LT*1000 - first_fib_ttl_ms; % Timestamp of TTL in ms relative to first ttl
cortex_pos_sliced =  cortex_pos(first_ttl_sample_ind_120hz:end, :);



first_ttl_sample_ind_100hz = round(first_fib_ttl_ms*ciholas_fs/1e3); % Sample index at 120hz corresponding to first ttl 
first_ttl_sample_ind_100hz = first_ttl_sample_ind_100hz - round(1.05*ciholas_fs);
last_ttl_sample_ind_100hz = round(last_fib_ttl_ms*ciholas_fs/1e3); % Sample index at 120hz corresponding to first ttl 
fib_ttl_ms = LT*1000 - first_fib_ttl_ms; % Timestamp of TTL in ms relative to first ttl
ciholas_pos_sliced =  extracted_cdp.tag_data_filt{1,1}(ciholas_first_sample_ind:end, 3:5);
ciholas_pos_120hz  = resample(ciholas_pos_sliced, 120, 100);
figure;
tiledlayout(3,1);
nexttile;
plot(cortex_pos_sliced(1:100000,1));
title('Cortex pos')
nexttile;
plot(ciholas_pos_120hz(1:100000,1));
title('Ciholas pos')
sgtitle('Cortex tracking seems to be missing a lot of data points. Perhaps we need add more markers to the cap')






