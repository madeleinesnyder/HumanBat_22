% Calibration done with test data from 220103
cortex = load('processed/220103/b149f/cortex/220103_14592_tracking_1_track.mat');
ciholas = load('processed/220103/b149f/ciholas/extracted_220103_cdp_1.mat');

ciholas_xpos_120hz = resample(ciholas.tag_data_filt{1}(:,3),120,100);
ciholas_ypos_120hz = resample(ciholas.tag_data_filt{1}(:,4),120,100);
ciholas_zpos_120hz = resample(ciholas.tag_data_filt{1}(:,5),120,100);

% ----------------------------------------------------------------
% Select start and end indices for segment of data to use for
% calibration:
ciholas_frame_s = round(20000); % Start of good calibration data with real movement trajectories
ciholas_frame_e = round(27000);
global_start_usec = ciholas.tag_data_filt{1}(ciholas_frame_s,8)*1e6;
global_end_usec = ciholas.tag_data_filt{1}(ciholas_frame_e,8)*1e6;
% ----------------------------------------------------------------
ciholas_frame_s_120hz = round(ciholas_frame_s*1.2); % Start of good calibration data with real movement trajectories
ciholas_frame_e_120hz = round(ciholas_frame_e*1.2);

% Slice ciholas position to window of interest.
ciholas_xpos_120hz = ciholas_xpos_120hz(ciholas_frame_s_120hz:ciholas_frame_e_120hz);
ciholas_ypos_120hz = ciholas_ypos_120hz(ciholas_frame_s_120hz:ciholas_frame_e_120hz);
ciholas_zpos_120hz = ciholas_zpos_120hz(ciholas_frame_s_120hz:ciholas_frame_e_120hz);

% Find corresponding index in cortex data for start and end as chosen in
% ciholas data
% Includes 1.05 second offset between cortex and ciholas
[minValue,start_index_cortex] = min(abs(cortex.global_sample_ts_usec-1.05*1e6 - global_start_usec));
[minValue,end_index_cortex] = min(abs(cortex.global_sample_ts_usec-1.05*1e6 - global_end_usec));
%cortex_xpos = HumanBat_interpolate_nans(cortex.avgMarkerPos(:,1));
%cortex_ypos = HumanBat_interpolate_nans(cortex.avgMarkerPos(:,2));
cortex_xpos = cortex.avgMarkerPos(start_index_cortex:end_index_cortex,1);
cortex_ypos = cortex.avgMarkerPos(start_index_cortex:end_index_cortex,2);
cortex_zpos = cortex.avgMarkerPos(start_index_cortex:end_index_cortex,3);

% Plot cortex positional data (120hz)
figure;
tiledlayout(4,1);
nexttile;
scatter(cortex_xpos, cortex_ypos);
nexttile;
plot(cortex_xpos);
nexttile;
plot(cortex_ypos);
nexttile;
plot(cortex_zpos);
sgtitle('Cortex trajectories 120hz')

% Plot 
figure;
tiledlayout(3,1);
nexttile;
scatter(ciholas_xpos_120hz, ciholas_ypos_120hz);
nexttile;
plot(ciholas_xpos_120hz);
nexttile;
plot(ciholas_ypos_120hz);
sgtitle('Ciholas trajectories 120hz')

% Ignore sample points with missing cortex tracking
cortex_good_samples = ~bitor(isnan(cortex_xpos), isnan(cortex_ypos));
cortex_good_samples = cortex_good_samples(1:min(length(ciholas_xpos_120hz), length(cortex_good_samples)));

% Only look at sample points with both ciholas and cortex tracking data
ciholas_x_samples = ciholas_xpos_120hz(cortex_good_samples);
ciholas_y_samples = ciholas_ypos_120hz(cortex_good_samples);
ciholas_z_samples = ciholas_zpos_120hz(cortex_good_samples);
cortex_x_samples = cortex_xpos(cortex_good_samples);
cortex_y_samples = cortex_ypos(cortex_good_samples);
cortex_z_samples = cortex_zpos(cortex_good_samples);

% Plot Sample points shared between cortex and ciholas
figure;
tiledlayout(6,1); 
nexttile;
plot(ciholas_x_samples);
title('ciholas x')
nexttile;
plot(cortex_x_samples);
title('cortex x')
nexttile;
plot(ciholas_y_samples);
title('ciholas y')
nexttile;
plot(cortex_y_samples);
title('cortex y')
nexttile;
plot(ciholas_z_samples);
title('ciholas z')
nexttile;
plot(cortex_z_samples);
title('cortex z')
sgtitle('Shared Non-NaN sample points')

% Get delay beteween cortex and ciholas. Should be 0!
x_delay = finddelay(cortex_x_samples, ciholas_x_samples);
y_delay = finddelay(cortex_y_samples, ciholas_y_samples);
z_delay = finddelay(cortex_z_samples, ciholas_z_samples);

fprintf('Should be close to 0: X delay: %d, Y delay: %d, Z delay: %d\n', x_delay, y_delay,z_delay);

% Get sample wise ratio between cortex and ciholas position
% Remove outliers.
cortex_ciholas_ratio = [cortex_x_samples./ciholas_x_samples cortex_y_samples./ciholas_y_samples cortex_z_samples./ciholas_z_samples];
cortex_ciholas_ratio = rmoutliers(cortex_ciholas_ratio);
x_cortex_ciholas_ratio = cortex_ciholas_ratio(:,1);
y_cortex_ciholas_ratio = cortex_ciholas_ratio(:,2);
z_cortex_ciholas_ratio = cortex_ciholas_ratio(:,3);

% Plot histogram of cortex/ciholas ratios. Should be close to normal
ratio_bins = linspace(0,1500,30);
figure;
tiledlayout(2,3);
nexttile;
histogram(x_cortex_ciholas_ratio, ratio_bins);
xlabel('ratio')
title('Cortex/Ciholas X pos ratio')
nexttile;
histogram(y_cortex_ciholas_ratio, ratio_bins);
xlabel('ratio')
title('Cortex/Ciholas Y pos ratio')
nexttile;
histogram(z_cortex_ciholas_ratio, ratio_bins);
xlabel('ratio')
title('Cortex/Ciholas Z pos ratio')
nexttile;
plot(x_cortex_ciholas_ratio);
xlabel('sample')
ylabel('ratio')
nexttile;
plot(y_cortex_ciholas_ratio);
xlabel('sample')
ylabel('ratio')
nexttile;
plot(z_cortex_ciholas_ratio);
xlabel('sample')
ylabel('ratio')

% Define conversion factor between ciholas and cortex as median of
% cortex/ciholas ratio
x_cortex_ciholas_factor = median(x_cortex_ciholas_ratio);
y_cortex_ciholas_factor = median(y_cortex_ciholas_ratio);
z_cortex_ciholas_factor = median(z_cortex_ciholas_ratio);

% Plot converted ciholas and cortex (should be on the same scale).
figure;
tiledlayout(3,1);
nexttile;
plot(ciholas_x_samples*x_cortex_ciholas_factor);
hold on;
plot(cortex_x_samples);
title('x')
nexttile;
plot(ciholas_y_samples*y_cortex_ciholas_factor);
hold on;
plot(cortex_y_samples);
title('y')
nexttile;
plot(ciholas_z_samples*z_cortex_ciholas_factor);
hold on;
plot(cortex_z_samples);
title('z')
sgtitle('Scaled cortex and ciholas positions')

% Save calibration for future use
ciholas2cortex = [x_cortex_ciholas_factor y_cortex_ciholas_factor z_cortex_ciholas_factor];
readme = 'Transformation factors based on data from 220103';
save('ciholas2cortex_scaling_factors.mat', 'ciholas2cortex', 'readme');