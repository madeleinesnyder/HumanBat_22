% This script is run via the batch file  hippo_extract_csc_batch_Michael.m .
% It extracts the following CSC-related data, and saves them into a mat-file
% (saves the data for each of the 3 sessions SEPARATELY!!):
% (1) EEG data: Extracts and saves the full-bandwidth (although Decimated) EEG data.
% (2) Theta-to-Delta Ratio: Computed in windows, with time-shifts that are smaller than the window size.
% (3) PSD (Power Spectral Density) in the same windows used to compute the Theta-to-Delta Ratio.
% (4) VELOCITY of the animal: Computed as 50'th or 90'th prctile (for the velocity)in the same windows used to compute the Theta-to-Delta Ratio.
% (5) Delta/Theta Oscillation Peak Times: Compute the times of the peak (findextrema.m) in the filtered EEG (filtered for Delta+Theta).
% (6) RIPPLES: Saves the Ripple Time, and the Ripple Shape.
% (7) Spikes fileNAMES: filenames of Acceptable single units ASSOCIATED with each CSC file.
%


% ======= Parameters: ===================================================
% Parameters are defined in the batch-script:  hippo_extract_csc_batch_Michael.m
% =======================================================================

%-----------------
% Michael Yartsev
%-----------------


disp(['Extracting CSC data, writing to:  ', filename_out]);


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
% --------------------------------------------------------------





% =======================================================================
% ======= Extract the SPIKE filenames associated with these CSC data (I do NOT require the units to be ACTIVE): =======
% =======================================================================

%load D:\Michael\Data\Expdata_Processed\Combined_data\data_inclusion_list_of_neurons.mat ; % Load the inclusion-list of Spike files

filename_associated_SPIKE_files = []; % Initialize
for ii_unit = 1:length(inclusion_list.file_list_allunits), % Loop over Spike filenames
    if ( strcmp( inclusion_list.file_list_associated_VT_files_allunits( ii_unit ), filename_associated_VT_file ) & ... % If this unit has the same VT file as the current CSC data (i.e. this is the same day)...
            inclusion_list.is_unit_is_pyramidal_is_active_allunits(ii_unit, 1) == 1 & ... % ...and this is a unit...
            inclusion_list.is_unit_is_pyramidal_is_active_allunits(ii_unit, 2) == 1 ), % ... and this is a pyramidal unit (but I do NOT require the unit to be ACTIVE)
        filename_associated_SPIKE_files = [ filename_associated_SPIKE_files , ...
            inclusion_list.file_list_allunits( ii_unit ) ] ; % Add this Spike file to the list
    end
end



% =======================================================================
% ======== Extract EEG data and compute the Theta-To-Delta Ratio + the Animal's Velocity: ============
% =======================================================================

clear  CSC ; % Clear this large and unnecessary variable

FieldSelection = [1 1 1 1 1] ; % Will read all the variables from the file
ExtractionMode = 1 ; % Will read all the file (sessions will be extracted later on)
ExtractionModeArray = [] ; % Will read all the file (sessions will be extracted later on)

[Timestamps, ChanNum, SampleFrequency, NumValidSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC( filename_CSC_EEG_in, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;

CSC_SamplePeriod_microsec = round( 10^6/SampleFrequency(1) );
CSC_Sampling_Rate_Hz = 10^6 / CSC_SamplePeriod_microsec ;
Samples_reshaped = zeros(1, prod(size(Samples)) );
Timestamps_filledIn = zeros(1, prod(size(Samples)) );
for ii_DataBlock = 1:size(Samples,2), % Loop over the 512-point blocks of data
    idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
    Samples_reshaped( idx_data ) = Samples(:,ii_DataBlock)';
    Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
end
Samples = Samples_reshaped - mean(Samples_reshaped); % The data, with Mean Removed
Timestamps = Timestamps_filledIn; % Timestamps in microsec

clear  Timestamps_filledIn  Samples_reshaped  ; % Clear some large and unnecessary variables, to avoid "OUT OF MEMORY" problems

% Resample (decimate) the data to a lower sampling rate:

CSC_SamplePeriod_microsec = CSC_SamplePeriod_microsec * decimate_factor ; % Re-define sampling rate/sampling period
CSC_Sampling_Rate_Hz = CSC_Sampling_Rate_Hz / decimate_factor ; % Re-define sampling rate/sampling period
Samples_decimated = decimate( Samples, decimate_factor ); % Decimate
Timestamps_decimated = (0:length(Samples_decimated)-1)*CSC_SamplePeriod_microsec + mean(Timestamps(1:2)); % Decimate; 1'st "decimated timestamp" was taken as the mean of the firat 2 original timestamps
% % % Samples_decimated([1 end]) = []; % Discard the end-points (decimation may be problematic there)
% % % Timestamps_decimated([1 end]) = []; % Discard the end-points (decimation may be problematic there)


% Convert the data ("Samples") to Microvolts and filter out (band-stop) the 50-Hz, 100-Hz & 150-Hz noises:

Samples_microvolts = Samples_decimated * str2num( NlxHeader{15}(13:end) ) * 10^6 ; % Convert to Microvolts
filter_Wn_Hz = [0.5 49 51 99 101 473] ; % For filtering OUT the 50Hz and 100Hz (the second harmonic) noise, and also the bands < 0.5 Hz and > 473 Hz
filter_N_order = 2000 ;  % For filtering OUT the 50Hz noise, and also the bands < 0.5 Hz and > 475 Hz
filter_filtertype = 'DC-0' ; % makes the first band (starting at 0 Hz) a stopband
filter_Wn = filter_Wn_Hz / ( CSC_Sampling_Rate_Hz / 2); % Normalize Wn to half the sampling rate
win = fir1( filter_N_order, filter_Wn, filter_filtertype ); % The multiband filtering window
%win_50Hz_bandstop = win_50Hz_bandstop / sum( win_50Hz_bandstop ); % NO NEED to normalize the fir1 filter by its sum!!!
Samples_microvolts_filtered50Hz = filtfilt( win, 1, Samples_microvolts ); % The filtering itself


% Compute Power-Spectral-Density using Welch's method (pwelch.m) with a Hamming window of
% N points ~2 sec (these numbers are based on: Csicsvari et al, J.Neurosci. 19, 274-287(1999)) --
% And then compute the theta-to-delta ratio ; also, compute the Velocity (50'th and 90'th prctile)) of the animal in the same segments (windows):

n_WelchWindows = floor( ( Timestamps_decimated(end) - Timestamps_decimated(1) ) / ThetaToDeltaRatio_TimeShiftWindow / 10^6 - ...
    ThetaToDeltaRatio_WindowLength / ThetaToDeltaRatio_TimeShiftWindow ) ; % Number of windows to loop over

ThetaToDeltaRatio_Ratio = zeros( n_WelchWindows, 1 ) + NaN ; % Initialize
ThetaToDeltaRatio_TimeOfWindowCenter = zeros( n_WelchWindows, 1 ) + NaN ;
ThetaToDeltaRatio_Raw_PSD = zeros( n_WelchWindows, ThetaToDeltaRatio_WelchWindowSize/2 + 1 ) + NaN ;
Velocity_Of_Animal_prctile50 = zeros( n_WelchWindows, 1 ) + NaN ;
Velocity_Of_Animal_prctile90 = zeros( n_WelchWindows, 1 ) + NaN ;
%Angular_Velocity_Of_Animal_mean = zeros( n_WelchWindows, 1 ) + NaN ;
Spike_Counts_of_individual_units = zeros( n_WelchWindows, length(filename_associated_SPIKE_files) ) + NaN ;

eval(['load ' filename_associated_VT_file]); % Load the VT file (for extracting the Velocity)

for ii_ThetaToDeltaRatio_window = 1:n_WelchWindows, % Loop over windows

    %ii_ThetaToDeltaRatio_window/n_WelchWindows
    ttt_segment = [ Timestamps_decimated(1) + (ii_ThetaToDeltaRatio_window - 1)*ThetaToDeltaRatio_TimeShiftWindow*10^6 , ...
        Timestamps_decimated(1) + (ii_ThetaToDeltaRatio_window - 1)*ThetaToDeltaRatio_TimeShiftWindow*10^6 + ThetaToDeltaRatio_WindowLength*10^6 ];

    % Compute the theta-to-delta ratio:
    idx_segment = find( Timestamps_decimated >= ttt_segment(1) & Timestamps_decimated <= ttt_segment(2) ); % Segment (window) for computing the theta-to-delta ratio

    [psd_eeg, freqs_eeg] = pwelch( Samples_microvolts_filtered50Hz( idx_segment ), ...
        ThetaToDeltaRatio_WelchWindowSize, [], [], CSC_Sampling_Rate_Hz ); % PSD over the Segment

    idx_Theta_freqs = find( freqs_eeg >= ThetaToDeltaRatio_ThetaFreqsHz(1) & ...
        freqs_eeg <= ThetaToDeltaRatio_ThetaFreqsHz(2) ); % Frequencies in the Theta range (Csicsvari et al, J.Neurosci. 19, 274-287(1999))
    idx_Delta_freqs = find( freqs_eeg >= ThetaToDeltaRatio_DeltaFreqsHz(1) & ...
        freqs_eeg <= ThetaToDeltaRatio_DeltaFreqsHz(2) ); % Frequencies in the Delta range (Csicsvari et al, J.Neurosci. 19, 274-287(1999))

    ThetaToDeltaRatio_Ratio( ii_ThetaToDeltaRatio_window ) = ...
        mean( psd_eeg( idx_Theta_freqs ) ) / mean( psd_eeg( idx_Delta_freqs ) ); % Ratio of the Mean PSD's in the 2 frequency ranges

    ThetaToDeltaRatio_TimeOfWindowCenter( ii_ThetaToDeltaRatio_window ) = ...
        Timestamps_decimated(1) + (ii_ThetaToDeltaRatio_window - 1)*ThetaToDeltaRatio_TimeShiftWindow*10^6 + ThetaToDeltaRatio_WindowLength*10^6 / 2 ; % Time of window-center (segment center) over which I computed the theta-to-delta ratio

    ThetaToDeltaRatio_Raw_PSD( ii_ThetaToDeltaRatio_window, : ) = psd_eeg(:)' ; % Save the raw PSD

    % Compute the animal's velocity in the same time segment (ONLY for the behavioral session):
    if ( ~isempty( VT ) ), % If this behavioral session was run at all
        if ( ttt_segment(1) >= VT.Timestamps(1) & ttt_segment(2) <= VT.Timestamps(end) ), % If this segment is inside the Behavioral Session (i.e. there IS velocity)
            idx_segment_VT = find( VT.Timestamps >= ttt_segment(1) & VT.Timestamps <= ttt_segment(2) ); % VT Segment for exracting the Velocity 
            Velocity_Of_Animal_prctile50( ii_ThetaToDeltaRatio_window ) = ... ; % 50'th percentile of the linear Velocity
                prctile( sqrt( sum(VT.velocity(idx_segment_VT,:)'.^2) ), 50 );
            Velocity_Of_Animal_prctile90( ii_ThetaToDeltaRatio_window ) = ... ; % 90'th percentile of the linear Velocity
                prctile( sqrt( sum(VT.velocity(idx_segment_VT,:)'.^2) ), 90 );
        end
    end

    % Compute the Number-of-Spikes within the same segment -- separately for each of the units:
    for ii_unit = 1:length(filename_associated_SPIKE_files), % Loop over units
        FieldSelection = [1 0 0 0 0] ; % Will read ONLY the Timestamps
        ExtractionMode = 1 ; % Mode 1 = "Extract All"
        ExtractionModeArray = [] ; % Will read all the data
        Timestamps_SPIKES = Nlx2MatSpike( filename_associated_SPIKE_files{ii_unit}, FieldSelection, 0, ExtractionMode, ExtractionModeArray ) ;
        idx_segment_SPIKES = find( Timestamps_SPIKES >= ttt_segment(1) & Timestamps_SPIKES <= ttt_segment(2) ); % Time Segment for Counting Spikes
        Spike_Counts_of_individual_units( ii_ThetaToDeltaRatio_window, ii_unit ) = length( idx_segment_SPIKES ); % Spike Count
    end

end % End "Loop over windows"



% Determine segments (windows) that include motion artefact (have very high power):

% --------
percent_for_trimmed_mean = 80 ; % Percent for the calculation of the trimmed mean below
% --------
mean_power_in_window = mean( ThetaToDeltaRatio_Raw_PSD' ) ; % Mean power in each of the windows
power_in_window_threshold = trimmean( mean_power_in_window, percent_for_trimmed_mean ) * ...
    ThetaToDeltaRatio_Mean_std_factor_for_detecting_artifacts ; % I used the TRIMMED mean trimmean() to avoid biasing the mean by the artifacts (outlier points) themseves...
ThetaToDeltaRatio_IS_OUTLIER_CHEWING_ARTIFACT = ( mean_power_in_window > ...
    power_in_window_threshold )'; % In the *PREVIOUS* line, I multiplied the mean power by the parameter "...Mean_std_factor...", to determine if the power in a segment (window) is a huge chewing-artifact

% Define/Save some more variables:

CSC_EEG_DECIMATED_SamplePeriod_microsec = CSC_SamplePeriod_microsec ;
CSC_EEG_DECIMATED_Sampling_Rate_Hz = CSC_Sampling_Rate_Hz;
Samples_microvolts_filtered50Hz_DECIMATED = Samples_microvolts_filtered50Hz ;
Timestamps_FirstTimestampOfEEGdata = Timestamps_decimated(1) ;
freqs_PSD_of_eeg = freqs_eeg ;




% =======================================================================
% ======= Find Delta/Theta Oscillation Peak Times (findextrema.m) =======
% =======================================================================

% Filter the EEG for Delta+Theta frequencies:

Wn = Delta_Plus_Theta__Filter_FreqsHz / ( CSC_Sampling_Rate_Hz / 2); % Passband, Hz (Normalized Wn to half the sampling rate)
win_ripples = fir1( Delta_Plus_Theta__Filter__N_order, Wn, 'bandpass' ); % Filter parameters are defined above (similar to Buzsaki)
%win_ripples = win_ripples / sum( win_ripples ); % No need to normalize the fir1 filter by its sum!!!
Samples_decimated_filt_DeltaTheta = filtfilt( win_ripples, 1, Samples_decimated ); % Use DECIMATED data

% Find the local maxima:

[Delta_Plus_Theta__idx_peaks,iLo,iCr] = findextrema(Samples_decimated_filt_DeltaTheta); 
Delta_Plus_Theta__peaks = Samples_decimated_filt_DeltaTheta(Delta_Plus_Theta__idx_peaks);
Delta_Plus_Theta__peaks_Timestamps = Timestamps_decimated( Delta_Plus_Theta__idx_peaks ); % Timestamps of the peaks
%
% figure ; hold on;
% plot( Timestamps_decimated, Samples_decimated_filt_DeltaTheta );
% plot( Delta_Plus_Theta__peaks_Timestamps, Delta_Plus_Theta__peaks, 'r.' );


% =======================================================================
% ======= Extract RIPPLE events: ==================
% =======================================================================

clear  Timestamps  Timestamps_filledIn  Timestamps_decimated  Samples  Samples_reshaped  Samples_decimated  Samples_microvolts  Samples_microvolts_filtered50Hz  Samples_decimated_filt_DeltaTheta ; % Clear some large and unnecessary variables, to avoid "OUT OF MEMORY" problems

% Extract the Samples/Timestamps data AGAIN, this time WITHOUT decimation:

FieldSelection = [1 1 1 1 1] ; % Will read all the variables from the file
ExtractionMode = 1 ; % Will read all the file (sessions will be extracted later on)
ExtractionModeArray = [] ; % Will read all the file (sessions will be extracted later on)

[Timestamps, ChanNum, SampleFrequency, NumValidSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC( filename_CSC_EEG_in, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;

CSC_SamplePeriod_microsec = round( 10^6/SampleFrequency(1) );
CSC_Sampling_Rate_Hz = 10^6 / CSC_SamplePeriod_microsec ;
Samples_reshaped = zeros(1, prod(size(Samples)) );
Timestamps_filledIn = zeros(1, prod(size(Samples)) );
for ii_DataBlock = 1:size(Samples,2), % Loop over the 512-point blocks of data
    idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
    Samples_reshaped( idx_data ) = Samples(:,ii_DataBlock)';
    Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
end
Samples = Samples_reshaped - mean(Samples_reshaped); % The data, with Mean Removed
Timestamps = Timestamps_filledIn; % Timestamps in microsec

% Save some more variables:

CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec = CSC_SamplePeriod_microsec ;

% Filter the EEG for Ripple frequencies:

clear  Timestamps_filledIn  Samples_reshaped ;  % Clear some large and unnecessary variables, to avoid "OUT OF MEMORY" problems

Wn = Ripples__Filter_FreqsHz / ( CSC_Sampling_Rate_Hz / 2); % Passband, Hz (Normalized Wn to half the sampling rate)
win_ripples = fir1( Ripples__Filter_N_order, Wn, 'bandpass' ); % Filter parameters are defined above (similar to Buzsaki)
%win_ripples = win_ripples / sum( win_ripples ); % No need to normalize the fir1 filter by its sum!!!
%%% Samples_filt_ripples = filtfilt( win_ripples, 1, Samples ); % Use NON-DECIMATED data
Samples_filt_ripples = []; % Initialize -- because of "Out Of Memory" problems, I will do the filtfilt() in Blocks
idx_block = 0 ; % Initialize -- to cover the (unlikely) situation that length(Samples) < n_points_block_size_filtfilt
n_points_block_size_filtfilt = 2*10^6; % I just chose a block size that is small enough NOT to cause "Out Of Memory" problems
for ii_block = 1 : floor( length(Samples) / n_points_block_size_filtfilt ),
    idx_block = (1:n_points_block_size_filtfilt)+(ii_block-1)*n_points_block_size_filtfilt ;
    Samples_filt_ripples = [ Samples_filt_ripples , ...
        filtfilt( win_ripples, 1, Samples( idx_block ) ) ];
end % The "problematic discontinuities" in Samples_filt_ripples - on the borders of 2 blocks - are SMALL in size
Samples_filt_ripples = [ Samples_filt_ripples , ... % Add the remainder of the data after the last block
    filtfilt( win_ripples, 1, Samples( idx_block(end)+1 : end ) ) ];

% Determine the samples which belong to Sleep Sessions (all the sleep sessions taken together!!!):

idx_sleep_sessions = []; % Initialize
if ( ~isempty( times_microsec_pre_sleep_session_1 ) ), % If this session exists
    idx_sleep_sessions = [idx_sleep_sessions, ...
        find( Timestamps >= times_microsec_pre_sleep_session_1(1) & ...
        Timestamps <= times_microsec_pre_sleep_session_1(2) ) ];
end
if ( ~isempty( times_microsec_post_sleep_session_3 ) ), % If this session exists
    idx_sleep_sessions = [idx_sleep_sessions, ...
        find( Timestamps >= times_microsec_post_sleep_session_3(1) & ...
        Timestamps <= times_microsec_post_sleep_session_3(2) ) ];
end

% Extract ONLY Sleep-Session data + Convert the Ripples "Samples" to MICROVOLTS, and then --
% Identify Ripples = EEG epochs that crossed the POWER threshold ; POWER = Abs( Hilbert ) ;
% Compute the Threshold ONLY over the SLEEP SESSIONS (all the sleep sessions taken together!!!):

Samples_filt_ripples = Samples_filt_ripples( idx_sleep_sessions ); % Cut only data belonging to Sleep Sessions!
Samples_filt_ripples = Samples_filt_ripples * str2num( NlxHeader{15}(13:end) ) * 10^6 ; % *** CONVERT TO MICROVOLTS ***
Timestamps = Timestamps( idx_sleep_sessions );

clear  Timestamps_filledIn  idx_sleep_sessions  Samples  Samples_reshaped  ; % Clear some large and unnecessary variables, to avoid "OUT OF MEMORY" problems

Samples_filt_ripples_POWER = []; % Initialize -- because of "Out Of Memory" problems, I will do the Hilbert transform in Blocks
idx_block = 0 ; % Initialize -- to cover the (unlikely) situation that length(Samples_filt_ripples) < n_points_block_size_hilbert
n_points_block_size_hilbert = 2*10^6; % I just chose a block size that is small enough NOT to cause "Out Of Memory" problems
for ii_block = 1 : floor( length(Samples_filt_ripples) / n_points_block_size_hilbert ),
    idx_block = (1:n_points_block_size_hilbert)+(ii_block-1)*n_points_block_size_hilbert ;
    Samples_filt_ripples_POWER = [ Samples_filt_ripples_POWER , ...
        abs( hilbert( Samples_filt_ripples( idx_block ) ) ) ];
end % The "problematic discontinuities" in Samples_filt_ripples_POWER - on the borders of 2 blocks - are SMALL in size
Samples_filt_ripples_POWER = [ Samples_filt_ripples_POWER , ... % Add the remainder of the data after the last block
    abs( hilbert( Samples_filt_ripples( idx_block(end)+1 : end ) ) ) ];

Ripples__power_threshold = mean( Samples_filt_ripples_POWER ) + ... % Determine the Threshold
    std( Samples_filt_ripples_POWER ) * Ripples__Threshold_std_factor ;

idx_Samples_Ripples_above_threshold = find( Samples_filt_ripples_POWER >= Ripples__power_threshold );

Timestamps_above_threshold = Timestamps( idx_Samples_Ripples_above_threshold ) ;

% Merge sub-portions of a Ripple that are divided by a SMALL gap (smaller than a
% threshold that I defined) -- and then extract the Timestamp of the max-point (PEAK) of the Ripple (this is PEAK for the LFP with NEGATIVITY POINTING UPWARDS).
% Also, extract the Ripples raw-data itself (Samples-filtered-for-ripples, extracted around the max-point of the Ripple):

idx_start_of_ripple = idx_Samples_Ripples_above_threshold(1) ; % Initialize

for ii_samples_over_threshold = 2:length(idx_Samples_Ripples_above_threshold), % Loop over above-threshold points

    if ( Timestamps_above_threshold(ii_samples_over_threshold) - ...
            Timestamps_above_threshold(ii_samples_over_threshold-1) > ...
            Ripples__time_gap_merging_threshold_ms * 10^3 ), % If there is a LONG time-gap from the PREVIOUS above-threshold sample -- define a NEW pulse

        idx_start_of_ripple = ... % idx of the NEW start-of-ripple
            [ idx_start_of_ripple , idx_Samples_Ripples_above_threshold( ii_samples_over_threshold ) ];

    end

end

% Loop over start-of-ripples, and extract the: (1) Ripples Timestamp = TIME-OF-PEAK (this is PEAK for the LFP with NEGATIVITY POINTING UPWARDS);
% (2) RIPPLE SHAPE (around the peak):

n_points_window_find_peak = round( Ripples__time_window_for_finding_peak_of_ripple_ms / ...
    CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec * 10^3 ); % Search for the ripple-peak inside a window defined by this time-window ON EACH SIDE of the start-of-ripple

n_points_window_save_ripple_waveform = round( Ripples__time_window_for_saving_ripple_waveform_ms / ...
    CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec * 10^3 ); % Save the ripple-waveform defined by this time-window ON EACH SIDE of the ripples peak

Ripples__Timestamp_peak_of_ripples = zeros( length(idx_start_of_ripple)-2, 1 ) + NaN ; % Initialize
Ripples__waveform_ripples = zeros( length(idx_start_of_ripple)-2, 2 * n_points_window_save_ripple_waveform + 1 );

for ii_ripple = 2:length(idx_start_of_ripple)-1, % Loop over Ripples, discarding the first and last ripple (I do this to prevent edge-problems)
    ripple_segment = Samples_filt_ripples( idx_start_of_ripple( ii_ripple ) - n_points_window_find_peak : ...
        idx_start_of_ripple( ii_ripple ) + n_points_window_find_peak );
    [stam, idx_max] = max( ripple_segment );
    idx_peak_of_ripple = idx_start_of_ripple( ii_ripple ) - n_points_window_find_peak + idx_max - 1 ;
    Ripples__Timestamp_peak_of_ripples( ii_ripple-1 ) = Timestamps( idx_peak_of_ripple );
    Ripples__waveform_ripples( ii_ripple-1, : ) = ...
        Samples_filt_ripples( idx_peak_of_ripple - n_points_window_save_ripple_waveform : ...
        idx_peak_of_ripple + n_points_window_save_ripple_waveform );
end


% Plot Ripples and their extraction:
% figure ; hold on
% plot( Timestamps, Samples_filt_ripples );
% plot( Timestamps, Samples_filt_ripples_POWER, 'r' );
% plot( Timestamps( idx_Samples_Ripples_above_threshold ), 0, 'g.' );
% plot( Timestamps( idx_start_of_ripple ), 10, 'k.' );
% plot(  Ripples__Timestamp_peak_of_ripples, 20, 'm.' );






% =======================================================================
% ======= Cut (some of) the data into SESSIONS: =======
% =======================================================================

times_session_allsessions{1} = times_microsec_pre_sleep_session_1; % Collet the times of all 3 sessions, so I could Loop on the sessions
times_session_allsessions{2} = times_microsec_behav_session_2;
times_session_allsessions{3} = times_microsec_post_sleep_session_3;

clear *___BySession* ; % Clear the old variables

for ii_session = 1 : 3, % Loop over all sessions

    times_session = times_session_allsessions{ ii_session }; % Times of current session

    if ( ~isempty( times_session ) ), % If this session WAS run

        % Cutting decimated-rate data within this session (but first I reconstituted the Timestamps_decimated variable from the 1'st sample's timestamp):
        Timestamps_decimated = [ 0 : length(Samples_microvolts_filtered50Hz_DECIMATED)-1 ] * ...
            CSC_EEG_DECIMATED_SamplePeriod_microsec + Timestamps_FirstTimestampOfEEGdata ;
        idx_BySession_decimated_rate = find( Timestamps_decimated >= times_session(1) & ...
            Timestamps_decimated <= times_session(2) );
        Samples_microvolts_filtered50Hz_DECIMATED___BySession{ ii_session } = ...
            Samples_microvolts_filtered50Hz_DECIMATED( idx_BySession_decimated_rate ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS
        Timestamps_FirstTimestampOfEEGdata___BySession{ ii_session } = ...
            Timestamps_decimated( idx_BySession_decimated_rate(1) ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS

        % Cutting window-PSD data within this session:
        idx_BySession_windows_PSD = find( ThetaToDeltaRatio_TimeOfWindowCenter > times_session(1) + ThetaToDeltaRatio_WindowLength*10^6 / 2  & ...
            ThetaToDeltaRatio_TimeOfWindowCenter < times_session(2) - ThetaToDeltaRatio_WindowLength*10^6 / 2 );
        ThetaToDeltaRatio_Ratio___BySession{ ii_session } = ...
            ThetaToDeltaRatio_Ratio( idx_BySession_windows_PSD ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS
        ThetaToDeltaRatio_TimeOfWindowCenter___BySession{ ii_session } = ...
            ThetaToDeltaRatio_TimeOfWindowCenter( idx_BySession_windows_PSD ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS
        ThetaToDeltaRatio_Raw_PSD___BySession{ ii_session } = ...
            ThetaToDeltaRatio_Raw_PSD( idx_BySession_windows_PSD, : ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS
        ThetaToDeltaRatio_IS_OUTLIER_CHEWING_ARTIFACT___BySession{ ii_session } = ...
            ThetaToDeltaRatio_IS_OUTLIER_CHEWING_ARTIFACT( idx_BySession_windows_PSD ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS
        Velocity_Of_Animal_prctile50___BySession{ ii_session } = ...
            Velocity_Of_Animal_prctile50( idx_BySession_windows_PSD ); % COLLECTING **BEHAVIORAL** SESSIONS ONLY
        Velocity_Of_Animal_prctile90___BySession{ ii_session } = ...
            Velocity_Of_Animal_prctile90( idx_BySession_windows_PSD ); % COLLECTING **BEHAVIORAL** SESSIONS ONLY
        Spike_Counts_of_individual_units___BySession{ ii_session } = ...
            Spike_Counts_of_individual_units( idx_BySession_windows_PSD, : ); % COLLECTING **BEHAVIORAL** SESSIONS ONLY

        % Cutting Delta+Theta Peaks data within this session:
        idx_BySession_ThetaDelta = find( Delta_Plus_Theta__peaks_Timestamps >= times_session(1) & ...
            Delta_Plus_Theta__peaks_Timestamps <= times_session(2) );
        Delta_Plus_Theta__peaks_Timestamps___BySession{ ii_session } = ...
            Delta_Plus_Theta__peaks_Timestamps( idx_BySession_ThetaDelta ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS
        Delta_Plus_Theta__peaks___BySession{ ii_session } = ...
            Delta_Plus_Theta__peaks( idx_BySession_ThetaDelta ); % COLLECTING SLEEP & BEHAVIORAL SESSIONS

        % Cutting Ripples data within this session:
        idx_BySession_ripples = find( Ripples__Timestamp_peak_of_ripples >= times_session(1) & ...
            Ripples__Timestamp_peak_of_ripples <= times_session(2) );
        Ripples__Timestamp_peak_of_ripples___BySession{ ii_session } = ...
            Ripples__Timestamp_peak_of_ripples( idx_BySession_ripples ); % COLLECTING **SLEEP** SESSIONS ONLY
        Ripples__waveform_ripples___BySession{ ii_session } = ...
            Ripples__waveform_ripples( idx_BySession_ripples, : ); % COLLECTING **SLEEP** SESSIONS ONLY


    else  % If this session was NOT run

        Samples_microvolts_filtered50Hz_DECIMATED___BySession{ ii_session } = 'Session was NOT run' ;
        Timestamps_FirstTimestampOfEEGdata___BySession{ ii_session } = 'Session was NOT run' ;
        ThetaToDeltaRatio_Ratio___BySession{ ii_session } = 'Session was NOT run' ;
        ThetaToDeltaRatio_TimeOfWindowCenter___BySession{ ii_session } = 'Session was NOT run' ;
        ThetaToDeltaRatio_Raw_PSD___BySession{ ii_session } = 'Session was NOT run' ;
        ThetaToDeltaRatio_IS_OUTLIER_CHEWING_ARTIFACT___BySession{ ii_session } = 'Session was NOT run' ;
        Velocity_Of_Animal_prctile50___BySession{ ii_session } = 'Session was NOT run' ;
        Velocity_Of_Animal_prctile90___BySession{ ii_session } = 'Session was NOT run' ;
        Spike_Counts_of_individual_units___BySession{ ii_session } = 'Session was NOT run' ;
        Delta_Plus_Theta__peaks_Timestamps___BySession{ ii_session } = 'Session was NOT run' ;
        Delta_Plus_Theta__peaks___BySession{ ii_session } = 'Session was NOT run' ;
        Ripples__Timestamp_peak_of_ripples___BySession{ ii_session } = 'Session was NOT run' ;
        Ripples__waveform_ripples___BySession{ ii_session } = 'Session was NOT run' ;

    end % End "If this session was run at all"

end % End "Loop over sessions"


% Put "NaN" in certain variables for certian sessions
% (e.g. in Velocity-Sleep-Sessions, or Ripples-Behavioral-Sessions):

[ Velocity_Of_Animal_prctile50___BySession{ [1 3] } ] = deal( NaN ); % COLLECTING **BEHAVIORAL** SESSIONS ONLY
[ Velocity_Of_Animal_prctile90___BySession{ [1 3] } ] = deal( NaN ); % COLLECTING **BEHAVIORAL** SESSIONS ONLY

[ Ripples__Timestamp_peak_of_ripples___BySession{ [2] } ] = deal( NaN ); % COLLECTING **SLEEP** SESSIONS ONLY
[ Ripples__waveform_ripples___BySession{ [2] } ] = deal( NaN ); % COLLECTING **SLEEP** SESSIONS ONLY


[ Timestamps_start_of_vocalizations___BySession{ [1 3] } ] = deal( NaN ); % COLLECTING **BEHAVIORAL** SESSIONS ONLY

% Clear some of the largest unnecessary varialbes:

clear  idx_*  Timestamps_decimated  Samples_microvolts_filtered50Hz_DECIMATED  ThetaToDeltaRatio_Raw_PSD  Ripples__waveform_ripples  Timestamps_above_threshold  Samples_above_threshold  Samples_AM_smoothed_above_threshold ;





% =======================================================================
% ======= Save the MAT-file with all the important variables: =======
% =======================================================================

clear CSC ;

CSC.params.decimate_factor = decimate_factor ;
CSC.params.ThetaToDeltaRatio_WindowLength = ThetaToDeltaRatio_WindowLength ;
CSC.params.ThetaToDeltaRatio_TimeShiftWindow = ThetaToDeltaRatio_TimeShiftWindow ;
CSC.params.ThetaToDeltaRatio_WelchWindowSize = ThetaToDeltaRatio_WelchWindowSize ;
CSC.params.ThetaToDeltaRatio_ThetaFreqsHz = ThetaToDeltaRatio_ThetaFreqsHz ;
CSC.params.ThetaToDeltaRatio_DeltaFreqsHz = ThetaToDeltaRatio_DeltaFreqsHz ;
CSC.params.ThetaToDeltaRatio_Mean_std_factor_for_detecting_artifacts = ThetaToDeltaRatio_Mean_std_factor_for_detecting_artifacts ;
CSC.params.Delta_Plus_Theta__Filter_FreqsHz = Delta_Plus_Theta__Filter_FreqsHz ;
CSC.params.Delta_Plus_Theta__Filter__N_order = Delta_Plus_Theta__Filter__N_order ;
CSC.params.Ripples__Filter_FreqsHz = Ripples__Filter_FreqsHz ;
CSC.params.Ripples__Filter_N_order = Ripples__Filter_N_order ;
CSC.params.Ripples__Threshold_std_factor = Ripples__Threshold_std_factor ;
CSC.params.Ripples__time_gap_merging_threshold_ms = Ripples__time_gap_merging_threshold_ms ;
CSC.params.Ripples__time_window_for_finding_peak_of_ripple_ms = Ripples__time_window_for_finding_peak_of_ripple_ms ;
CSC.params.Ripples__time_window_for_saving_ripple_waveform_ms = Ripples__time_window_for_saving_ripple_waveform_ms ;
CSC.params.filename_CSC_EEG_in = filename_CSC_EEG_in ;
CSC.params.filename_associated_VT_file = filename_associated_VT_file ;
CSC.params.times_microsec_pre_sleep_session_1 = times_microsec_pre_sleep_session_1 ;
CSC.params.times_microsec_behav_session_2 = times_microsec_behav_session_2 ;
CSC.params.times_microsec_post_sleep_session_3 = times_microsec_post_sleep_session_3 ;

CSC.data.Samples_microvolts_filtered50Hz_DECIMATED___BySession = Samples_microvolts_filtered50Hz_DECIMATED___BySession ;
CSC.data.Timestamps_FirstTimestampOfEEGdata___BySession = Timestamps_FirstTimestampOfEEGdata___BySession ;
CSC.data.CSC_EEG_DECIMATED_SamplePeriod_microsec = CSC_EEG_DECIMATED_SamplePeriod_microsec ;
CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec = CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec ;
CSC.data.CSC_EEG_DECIMATED_Sampling_Rate_Hz = CSC_EEG_DECIMATED_Sampling_Rate_Hz;
CSC.data.ThetaToDeltaRatio_Ratio___BySession = ThetaToDeltaRatio_Ratio___BySession ;
CSC.data.ThetaToDeltaRatio_TimeOfWindowCenter___BySession = ThetaToDeltaRatio_TimeOfWindowCenter___BySession ;
CSC.data.ThetaToDeltaRatio_Raw_PSD___BySession = ThetaToDeltaRatio_Raw_PSD___BySession ;
CSC.data.ThetaToDeltaRatio_IS_OUTLIER_CHEWING_ARTIFACT___BySession = ThetaToDeltaRatio_IS_OUTLIER_CHEWING_ARTIFACT___BySession ;
CSC.data.freqs_PSD_of_eeg = freqs_PSD_of_eeg ;
CSC.data.Delta_Plus_Theta__peaks_Timestamps___BySession = Delta_Plus_Theta__peaks_Timestamps___BySession ;
CSC.data.Delta_Plus_Theta__peaks___BySession = Delta_Plus_Theta__peaks___BySession ;
CSC.data.Velocity_Of_Animal_prctile50___BySession = Velocity_Of_Animal_prctile50___BySession ;
CSC.data.Velocity_Of_Animal_prctile90___BySession = Velocity_Of_Animal_prctile90___BySession ;
CSC.data.Spike_Counts_of_individual_units___BySession = Spike_Counts_of_individual_units___BySession ;
CSC.data.Ripples__power_threshold = Ripples__power_threshold ;
CSC.data.Ripples__Timestamp_peak_of_ripples___BySession = Ripples__Timestamp_peak_of_ripples___BySession ;
CSC.data.Ripples__waveform_ripples___BySession = Ripples__waveform_ripples___BySession ;
CSC.data.tetrode_depth_nominal_microns = tetrode_depth_nominal_microns ;
CSC.data.filename_associated_SPIKE_files = filename_associated_SPIKE_files ;

eval( ['save ', filename_out, ' CSC'] );


% ================

cd D:\Michael\Matlab ;

% ======= THE END ===============

