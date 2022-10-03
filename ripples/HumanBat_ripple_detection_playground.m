for i = 1:4
    [Timestamps, SampleFrequency, c] = getRawCSCData(sprintf('CSC%d.ncs',i),1,10000000);
    csc(:,i) = c;
end

% ------------  General Parameters: ----------------------------
decimate_factor = 2 ; % Decimate the EEG by this factor prior to computing Welch's PSD
ThetaToDeltaRatio_WindowLength = 2.18 ; % Window size (sec) for computing the Theta-to-Delta Ratio -- MUST be LONGER than ThetaToDeltaRatio_WelchWindowSize (BTW, this is MUCH shorter than the 15-sec MINIMAL time of a "theta-epoch" mentioned in: Buzsakei et al, Neuroscience 116, 201-211(2003))
ThetaToDeltaRatio_TimeShiftWindow = 1 ; % Time shifts of Windows (sec) when computing the Theta-to-Delta Ratio
ThetaToDeltaRatio_WelchWindowSize = 2048 ; % Window size (samples) for power-spectral-density -- MUST be SHORTER than ThetaToDeltaRatio_WindowLength (indeed it is ~2 sec, similar to: Csicsvari et al, J.Neurosci. 19, 274-287(1999))
ThetaToDeltaRatio_ThetaFreqsHz = [4 8] ; % Theta freqs (Hz) used for computing the Theta-to-Delta Ratio (numbers are LOWER than: Csicsvari et al, J.Neurosci. 19, 274-287(1999))
ThetaToDeltaRatio_DeltaFreqsHz = [2 4] ;  % Delta freqs (Hz) used for computing the Theta-to-Delta Ratio (numbers are LOWER than: Csicsvari et al, J.Neurosci. 19, 274-287(1999))
ThetaToDeltaRatio_Mean_std_factor_for_detecting_artifacts = 10 ; % Will multiply the mean power by this factor to determine that the power in a segment (window) is a huge chewing-artifact
Delta_Plus_Theta__Filter_FreqsHz = [0.1 9] ; % Filtering Delta+Theta frequencies (for finding the peaks of the oscillations)
Delta_Plus_Theta__Filter__N_order = 2000 ; % Filtering Delta+Theta frequencies (for finding the peaks of the oscillations)
Ripples__Filter_FreqsHz = [80 160]; % Detecting Ripples: Frequencies for filtering the EEG for Ripples [numbers are LOWER than: Buzsakei et al, Neuroscience 116, 201-211(2003)]
Ripples__Filter_N_order = 30 ; % Detecting Ripples: Order of filter (filtering for Ripples) [numbers from: Buzsakei et al, Neuroscience 116, 201-211(2003)]
Ripples__Threshold_std_factor = 7 ; % A Ripples is defined when POWER of filtered-EEG > mean + std*(this factor) [numbers from: Buzsakei et al, Neuroscience 116, 201-211(2003)]
Ripples__time_gap_merging_threshold_ms = 25 ; % I will merge Ripples separated by a time-gap smaller than this number (ms)
Ripples__time_window_for_finding_peak_of_ripple_ms = 30 ; % Search for the ripple-peak inside a window defined by this time-window ON EACH SIDE of the start-of-ripple
Ripples__time_window_for_saving_ripple_waveform_ms = 200 ; % Save the ripple-waveform defined by this time-window ON EACH SIDE of the ripples peak
CSC_SamplePeriod_microsec = round( 10^6/SampleFrequency(1) );
CSC_Sampling_Rate_Hz = SampleFrequency(1)%10^6 / CSC_SamplePeriod_microsec ;
%Samples_reshaped = zeros(1, prod(size(Samples)) );
% --------------------------------------------------------------

Samples = csc(:,1);
Timestamps_filledIn = zeros(1, prod(size(Samples)) );
for ii_DataBlock = 1:size(Samples,2), % Loop over the 512-point blocks of data
    idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
    Samples_reshaped( idx_data ) = Samples(:,ii_DataBlock)';
    Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
end
Samples = Samples_reshaped - mean(Samples_reshaped); % The data, with Mean Removed
Timestamps = Timestamps_filledIn; % Timestamps in microsec

CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec = CSC_SamplePeriod_microsec ;

Wn = Ripples__Filter_FreqsHz / ( CSC_Sampling_Rate_Hz / 2); % Passband, Hz (Normalized Wn to half the sampling rate)
win_ripples = fir1( Ripples__Filter_N_order, Wn, 'bandpass' ); % Filter parameters are defined above (similar to Buzsaki)
%win_ripples = win_ripples / sum( win_ripples ); % No need to normalize the fir1 filter by its sum!!!
%%% Samples_filt_ripples = filtfilt( win_ripples, 1, Samples ); % Use NON-DECIMATED data
Samples_filt_ripples = []; % Initialize -- because of "Out Of Memory" problems, I will do the filtfilt() in Blocks
idx_block = 0 ; % Initialize -- to cover the (unlikely) situation that length(Samples) < n_points_block_size_filtfilt
n_points_block_size_filtfilt = 2*10^6; % I just chose a block size that is small enough NOT to cause "Out Of Memory" problems
Samples_filt_ripples = bandpass(Samples, Ripples__Filter_FreqsHz, CSC_Sampling_Rate_Hz);
%Samples_filt_ripples = filtfilt(win_ripples, 1, Samples);

Samples_filt_ripples_POWER = abs(hilbert(Samples_filt_ripples));

Ripples__power_threshold = mean( Samples_filt_ripples_POWER ) + ... % Determine the Threshold
    std( Samples_filt_ripples_POWER ) * Ripples__Threshold_std_factor ;

idx_Samples_Ripples_above_threshold = find( Samples_filt_ripples_POWER >= Ripples__power_threshold );
Timestamps_above_threshold = Timestamps( idx_Samples_Ripples_above_threshold ) ;

figure;
tiledlayout(3,1);
nexttile
plot(Samples)
nexttile
plot(Samples_filt_ripples)
nexttile
plot(Samples_filt_ripples_POWER)

idx_start_of_ripple = idx_Samples_Ripples_above_threshold(1) ; % Initialize

for ii_samples_over_threshold = 2:length(idx_Samples_Ripples_above_threshold), % Loop over above-threshold points

    if ( Timestamps_above_threshold(ii_samples_over_threshold) - ...
            Timestamps_above_threshold(ii_samples_over_threshold-1) > ...
            Ripples__time_gap_merging_threshold_ms * 10^3 ), % If there is a LONG time-gap from the PREVIOUS above-threshold sample -- define a NEW pulse

        idx_start_of_ripple = ... % idx of the NEW start-of-ripple
            [ idx_start_of_ripple , idx_Samples_Ripples_above_threshold( ii_samples_over_threshold ) ];

    end

end
figure;
tiledlayout(3,4)
for i = 1:length(idx_start_of_ripple)
    ripple_idx = idx_start_of_ripple(i);
    s = ripple_idx - 5*30*32000/1000;
    e = ripple_idx + 10*30*32000/1000;
    xline(5*30*32000/1000,'k');
    axs(i) = nexttile
    plot(Samples(s:e));
    xlim([0 16000])
    xticks([0 8000 16000])
    xticklabels({'0','250','500'})
    xlabel('time (ms)')
    ylabel('CSC')
end
sgtitle('Detected ripples CSC (line denotes detected ripple start)')
linkaxes(axs, 'xy');

sampling_rate_Hz = 32000;
PSD_Ripples_npoints = 1024*12 ; % npoints for computing the PSD of LFP around Ripples ( this is for NONDECIMATED LFP DATA)
window_length_cut_sec = PSD_Ripples_npoints / sampling_rate_Hz * 1.05 ;

[total_pxx, total_fs] = pwelch(Samples, PSD_Ripples_npoints, [], [], sampling_rate_Hz ); % PSD, Welch method

ripple_segments = [] % Concatenated CSC of segments of csc traces in a 1 sec time window around ripples
for i = 1:length(idx_start_of_ripple)
    ripple_idx = idx_start_of_ripple(i);
    s = ripple_idx - 500*32000/1000; % 500ms before
    e = ripple_idx + 500*32000/1000; % 500ms after
    ripple_segments = [ripple_segments Samples(s:e)]
end
[ripple_pxx, ripple_fs] = pwelch(ripple_segments, PSD_Ripples_npoints, [], [], sampling_rate_Hz ); % PSD, Welch method

figure;
plot(total_fs, total_pxx,'DisplayName','Overall LFP power spectrum', 'Color','k')
set(gca, 'YScale', 'log')
hold on
plot(ripple_fs, ripple_pxx,'DisplayName','Average ripple segments power spectrum', 'Color', 'red')
xlim([0 240])
xticks([0 60 120 180 240])
xlabel('Frequency (hz)');
ylabel('Power (log)');

figure;
tiledlayout(3,4)
for i = 1:length(idx_start_of_ripple)
    ripple_idx = idx_start_of_ripple(i);
    s = ripple_idx - 1*30*32000/1000;
    e = ripple_idx + 2*30*32000/1000;

    [pxx, fs] = pwelch(Samples(s:e), PSD_Ripples_npoints, [], [], sampling_rate_Hz ); % PSD, Welch method
    
    axs(i) = nexttile
    plot(fs,pxx)
    hold on
    plot(total_fs, total_pxx)
    set(gca, 'YScale', 'log')
    xlim([0 240])
    xticks([0 60 120 180 240])
    xlabel('Frequency (hz)');
    ylabel('Power (log)');
end
sgtitle('Detected ripples PSD')
linkaxes(axs, 'xy');




% Compute the PSD for ALL the data in the Sleep sessions:

