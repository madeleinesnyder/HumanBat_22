% Analyzing ripples in my recordings -- plots some figures -- and saves the data into a mat-file:
% D:\Michael\Data\Expdata_Processed\Combined_data\data_theta.mat .

%-----------------
% Michael Yartsev
%-----------------


clear all ; close all ; fclose all ; pack ;


% ======= Parameters (NO BATCH FILE IS USED HERE!!): ==================
% filename_in_list = { ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day12\CSC_extracted_bat6053_Day12_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day13\CSC_extracted_bat6053_Day13_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day14\CSC_extracted_bat6053_Day14_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day17\CSC_extracted_bat6053_Day17_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day21\CSC_extracted_bat6053_Day21_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day22\CSC_extracted_bat6053_Day22_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day23\CSC_extracted_bat6053_Day23_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day24\CSC_extracted_bat6053_Day24_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day25\CSC_extracted_bat6053_Day25_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day26\CSC_extracted_bat6053_Day26_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day29\CSC_extracted_bat6053_Day29_tt3', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day13\CSC_extracted_bat7545_Day13_tt4', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day20\CSC_extracted_bat7545_Day20_tt4', ...
%     'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day21\CSC_extracted_bat7545_Day21_tt4', ...
%     };

% Analyzing days for which we saw place cells
filename_in_list = { ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day13\CSC_extracted_bat6053_Day13_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day14\CSC_extracted_bat6053_Day14_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day17\CSC_extracted_bat6053_Day17_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day22\CSC_extracted_bat6053_Day22_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day23\CSC_extracted_bat6053_Day23_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day24\CSC_extracted_bat6053_Day24_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day25\CSC_extracted_bat6053_Day25_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day26\CSC_extracted_bat6053_Day26_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day29\CSC_extracted_bat6053_Day29_tt3_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day13\CSC_extracted_bat7545_Day13_tt4_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day16\CSC_extracted_bat7545_Day16_tt4_flying_bouts_removed', ...
    'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day20\CSC_extracted_bat7545_Day20_tt4_flying_bouts_removed', ...
    };

% ======================================================================


% ------------  General Parameters: ----------------------------
PSD_Ripples_npoints = 1024 ; % npoints for computing the PSD of LFP around Ripples ( this is for DECIMATED LFP DATA)
PSD_Ripples_freqs_50Hz_harmonics_to_interpolate_in_plots = [48 52.5 ; 98 102.5 ; 148 152 ; 198 202] ; % Freqs of 50 Hz noises (and harmonics) to interpolate when plotting
PSD_Ripples_xticks = 0 : 50 : 1000 ; % Xticks for plotting PSD of ripples
PSD_Ripples_xlim = [0 320]; % Xlim for plotting PSD of ripples
Sharpwave_extract_npoints = 1000 ; % npoints for extracting the sharpwave around ripples( this is for DECIMATED LFP DATA)
Sharpwave_filter_Wn_Hz = [0.1 20] ; % Filtering the LFP for Sharpwaves
Sharpwave_filter_N_order = 200 ; % Filtering the LFP for Sharpwaves
xlim_sharpwave = [-500 500] ; % Xlimits for plotting the average Sharpwave (ms)
percent_trimmean_for_computing_mean_Ripple_size = 98 ; % Percent for computing trimmean of Ripple size
window_compute_ripple_size_sharpwave_size__ms = [-20 20]; % Time window (around ripple peak) for computing ripple size and sharpwave size
bins_microns_for_computing_avg_population_ripple_size = [ -180 -40 ; -40 40 ; 40 180 ]; % Bins (microns) for computing population-average of Ripple size
window_size_psth_for_computing_around_ripple_ms = 1000 ; % Window size for computing PSTHs around ripples (half this window on each size)
npoints_smooth_psth = 3 ; % Number of points ( = number of milliseconds) for smoothing the ripples-triggered PSTH
xlim_psth = [-300 300]; % Xlimits for plotting the ripples-triggered PSTH (ms)
xticks_psth = -500 : 100 : 500 ; % Xticks of the ripples-triggered PSTH (ms)
xlim_psth_ZoomIn = [-100 100]; % Xlimits of the ripples-triggered PSTH (ms)
xticks_psth_ZoomIn = -100 : 20 : 100 ; % Xticks of the ripples-triggered PSTH (ms)
dir_save_figs = 'D:\Michael\Data\Expdata_Processed\Combined_data\2_D_Hippo\Ripples'; % Here I will save the figures
idx_days_bat_6053 = 1:16 ; % Days recorded from ONE of the bats
idx_days_bat_7545 = 17:27 ; % Days recorded from ONE of the bats
idx_days_first2weeks_use_for_depth_profile = [1:16 17:27]  ; % Use only the first 2 weeks of data for depth-profile plotting (the depth at later dates may be not reliable)
InterRipplesInterval__discard_intervals_above_this_number_ms = 15000 ; % Discard intervals above this interval (ms)
%InterRipplesInterval__bins_for_histogram_ms = 22 : 44 : 5000 ; % Bins for computing inter-ripple-interval histogram
InterRipplesInterval__bins_for_histogram_ms = 0 : 100 : 5000 ; % Bins for computing inter-ripple-interval histogram
InterRipplesInterval__xlim_ms = [-30 3000] ; % Xlimits when plotting inter-ripples-intervals
minimal_num_ripples = 20; % Mininal number of ripples to inlcude a day into some of the analysis
filename_mat_file_save_data = 'D:\Michael\Data\Expdata_Processed\Combined_data\data_ripples_hippo.mat' ; % I will save the data into this mat-file
% --------------------------------------------------------------




PSD_Around_Ripples_SUM__all = []; % Initialize
PSD_Sleep_ALL__all = []; 
freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all = [];
Ripples_waveforms_mean__all = [];
psth_triggered_around_ripples__all_units__sp2sec__all = [];
psth_triggered_around_ripples__all_units__sp2sec__smoothed__all = [];
psth_triggered_around_ripples__all_units__smoothed__norm__all = [];
psth_triggered_around_ripples__multi_units__sp2sec__all = [];
psth_triggered_around_ripples__multi__sp2sec__smoothed__all = [];
psth_triggered_around_ripples__multi__smoothed__norm__all = [];
Depth_of_tetrode_nominal__all = [];
Day_of_single_unit_all = [];
Ripple_peak_to_peak_microvolts__Individual_Ripples__all = [];
Ripple_peak_to_peak_microvolts__mean__all = [];
Ripple_peak_to_peak_microvolts__sem__all = [];
Ripple_DayNumber_for_each_Individual_Ripple__all = [];
Ripple_Tetrode_Depth_for_each_Individual_Ripple__all = [];
Sharpwave_peak_to_peak_microvolts__mean__all = [];
Sharpwave_peak_to_peak_microvolts__sem__all = [];
InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms__all = [] ;
InterRippleIntervals_Hist_normalized_to_percent_occurrence__all = [] ;
Ripple_Duration_ms__all = [] ;
Ripple_frequency_Computed_using_local_maxima_all = [];
Timestamps_Ripples__SavedSeparatelyForEachDay = {} ;
Number_of_Ripples = [];


        
for ii_recording_day = 1:length(filename_in_list), % Loop over recording-days

    
    h = figure ; 
    % Some WYSIWYG options:
    set(gcf,'DefaultAxesFontSize',10);
    set(gcf,'DefaultAxesFontName','helvetica');
    set(gcf,'PaperUnits','centimeters','PaperPosition',[2 3 29 10]);
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0.2 0.2 0 0]);
    whitebg('black')

    
        
    clear  CSC  ; % Clear the largest variables, to free up space in memory
    
    eval(['load ', filename_in_list{ii_recording_day}]); % Load CSC data
    
    disp(['Analyzing RIPPLES from day # ', num2str(ii_recording_day), ' -- ', filename_in_list{ii_recording_day}]);
    

           
        
         
    % ==============================================================
    % ======== Plot PSD of Ripples -- Computed for All the sleep session, and also just around the Ripples: ============
    % ==============================================================
    
    
    % Loop over Behavioral Sessions, and read Samples_decimated and Timestamps_decimated (Timestamps_decimated
    % reconstituted from the 1'st-sample-timestamp), as well as the Timestamps of Ripples, 
    % as well as the Ripple waveforms:

    Samples_LFP_Sleep = []; % Initialize
    Timestamps_LFP_Sleep = [];
    Timestamps_Ripples = [];
    Ripples_waveforms = [];
    for ii_session = [1 3], % Loop over Sleep sessions
        if ( length( CSC.data.Samples_microvolts_filtered50Hz_DECIMATED___BySession{ii_session} ) > 99.9 & ... % If this session was run at all
             length(CSC.data.Ripples__Timestamp_peak_of_ripples___BySession{ii_session}) >= 1 ), % & if there were Ripple events in this session
                Samples_LFP_Sleep = [ Samples_LFP_Sleep , ...
                    CSC.data.Samples_microvolts_filtered50Hz_DECIMATED___BySession{ii_session} ];
                Timestamps_LFP_Sleep = [ Timestamps_LFP_Sleep , ...
                    CSC.data.Timestamps_FirstTimestampOfEEGdata___BySession{ii_session} + ...
                    [ 0 : length( CSC.data.Samples_microvolts_filtered50Hz_DECIMATED___BySession{ii_session} )-1 ] * ...
                    CSC.data.CSC_EEG_DECIMATED_SamplePeriod_microsec ]; % Reconstituted the Timestamps from the 1'st-sample-timestamp
                Timestamps_Ripples = [ Timestamps_Ripples ; CSC.data.Ripples__Timestamp_peak_of_ripples___BySession{ii_session} ]; % Collecting Timestamps of Ripples
                Ripples_waveforms = [ Ripples_waveforms ; CSC.data.Ripples__waveform_ripples___BySession{ii_session} ]; % Collecting Waveforms of Ripples
         end % End "If this session was run at all & if Ripples events occurred"
    end % End "Loop over Behavioral sessions"

    % Read LFP data in a T-sec window around each Ripple (T = window_length_cut_sec), and compute the PSD from these windows:
    
    PSD_Around_Ripples_SUM = 0; % Initialize
    for ii_ripple = 1:length(Timestamps_Ripples), % Loop over Ripples
        sampling_rate_Hz = 10^6 / CSC.data.CSC_EEG_DECIMATED_SamplePeriod_microsec ; % Sampling rate of the LFP (DECIMATED)
        window_length_cut_sec = PSD_Ripples_npoints / sampling_rate_Hz * 1.05 ;
        idx_data_for_PSD = find( abs( Timestamps_LFP_Sleep - Timestamps_Ripples(ii_ripple) ) <= window_length_cut_sec/2 * 10^6 ); % Cut a "window" (chunk) of LFP around each Ripple
        if ( length( idx_data_for_PSD ) >= PSD_Ripples_npoints ), % If the ripple is not too close to a session's start/end (which will truncate the "window" that I'm cutting)
            [PSD_Around_Ripples, PSD_freqs_Ripples] = pwelch( Samples_LFP_Sleep( idx_data_for_PSD ) , ...
                PSD_Ripples_npoints, [], [], sampling_rate_Hz ); % PSD, Welch method
            PSD_Around_Ripples_SUM = PSD_Around_Ripples_SUM + PSD_Around_Ripples / length(Timestamps_Ripples) ; % Summing (averaging) the PSD's
        end
        %subplot(4,5,ii_ripple); hold on; plot( Timestamps_LFP_Sleep( idx_data_for_PSD ), Samples_LFP_Sleep( idx_data_for_PSD ) ); plot( Timestamps_Ripples(ii_ripple), 0, 'r.' );
    end % End Loop over Ripples"

    % Compute the PSD for ALL the data in the Sleep sessions:
    
    [PSD_Sleep_ALL, PSD_freqs_Ripples] = pwelch( Samples_LFP_Sleep , ...
        PSD_Ripples_npoints, [], [], sampling_rate_Hz ); % PSD, Welch method
    
    % Remove PSD_Ripples_freqs_50 Hz noises (filtered Notches, or Peaks) from the PSD:
    
    for ii_PSD_Ripples_freqs_50Hz_harmonics = 1:size(PSD_Ripples_freqs_50Hz_harmonics_to_interpolate_in_plots,1), % Loop over PSD_Ripples_freqs_50 Hz harmonics
        idx_freqs_to_interpolate = find( PSD_freqs_Ripples >= PSD_Ripples_freqs_50Hz_harmonics_to_interpolate_in_plots( ii_PSD_Ripples_freqs_50Hz_harmonics, 1 ) & ...
            PSD_freqs_Ripples <= PSD_Ripples_freqs_50Hz_harmonics_to_interpolate_in_plots( ii_PSD_Ripples_freqs_50Hz_harmonics, 2 ) );
        PSD_Sleep_ALL( idx_freqs_to_interpolate ) = ...
            interp1( PSD_freqs_Ripples( [ idx_freqs_to_interpolate(1)-1  idx_freqs_to_interpolate(end)+1 ] ), ...
                PSD_Sleep_ALL( [ idx_freqs_to_interpolate(1)-1  idx_freqs_to_interpolate(end)+1 ] ), ...
                PSD_freqs_Ripples( idx_freqs_to_interpolate ), 'linear' ); % Linear interpolation using the 2 data-points that flank the frequency-range that is being cut
        PSD_Around_Ripples_SUM( idx_freqs_to_interpolate ) = ...
            interp1( PSD_freqs_Ripples( [ idx_freqs_to_interpolate(1)-1  idx_freqs_to_interpolate(end)+1 ] ), ...
                PSD_Around_Ripples_SUM( [ idx_freqs_to_interpolate(1)-1  idx_freqs_to_interpolate(end)+1 ] ), ...
                PSD_freqs_Ripples( idx_freqs_to_interpolate ), 'linear' ); % Linear interpolation using the 2 data-points that flank the frequency-range that is being cut
    end
    
    % Turn into row vectors:
    PSD_Sleep_ALL = PSD_Sleep_ALL(:)';
    PSD_Around_Ripples_SUM = PSD_Around_Ripples_SUM(:)';
    PSD_freqs_Ripples = PSD_freqs_Ripples(:)';
    
    
    % Compute the Difference between the two PSD's, and find freq( max of this difference ) = Freq of Ripple Peak:
    diff_of_PSDs_logscale = log10( PSD_Around_Ripples_SUM ) - log10( PSD_Sleep_ALL );
    [max_value, idx_max] = max( diff_of_PSDs_logscale );
    freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep = PSD_freqs_Ripples( idx_max );
    

    N_Ripples = length(Timestamps_Ripples); % Number of ripples we recorded today

    % Plot:
    subplot(2,3,[1,4]); hold on ;
        ylimits = log10( [ min(PSD_Sleep_ALL(PSD_Ripples_xlim(2)), PSD_Around_Ripples_SUM(PSD_Ripples_xlim(2)))   max([PSD_Sleep_ALL(:) ; PSD_Around_Ripples_SUM(:)]) ] );
        plot( PSD_freqs_Ripples, diff_of_PSDs_logscale, '-', 'color', 'c', 'linewidth', 0.5 ); 
        plot( freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep, max_value, 'r.', 'markersize', 8 ); % Marking the peak of hte Ripple
        plot( PSD_freqs_Ripples, log10( PSD_Sleep_ALL ), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5 );
        plot( PSD_freqs_Ripples, log10( PSD_Around_Ripples_SUM ), '-', 'color', 'g', 'linewidth', 1.5 );   
        title( { filename_in_list{ii_recording_day}(42:48), filename_in_list{ii_recording_day}(63:67)  ['Green:  PSD around RIPPLES ,     N ripples = ', num2str(length(Timestamps_Ripples))], ...
            'Gray: ALL the Sleep data ; Cyan = Diff', ...
            ['Peak Freq of Ripple (Hz) = ', num2str(freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep,4)]}, ...
            'interpreter', 'none' );
        xlabel( 'Frequency (Hz)' );
        ylabel('log_1_0 Power');
        set( gca, 'xlim', PSD_Ripples_xlim, 'ylim', ylimits, 'xtick', PSD_Ripples_xticks, 'box', 'off' );

    
        
    
    
    % =================================================================================
    % ======== Plot the average Ripple Shape + Ripple Power (Abs-Hilbert): ============
    % =================================================================================

    % Some useful data -- will use this for the Depth-profiles later on:
    
    Depth_of_tetrode_nominal = CSC.data.tetrode_depth_nominal_microns ; % Will use this for the Depth-profiles later on:

    % Plot:
    
        Ripples_waveforms_mean = mean( Ripples_waveforms );
        ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
            ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
            CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
            10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
        idx_for_computing_ripples_size_and_sharpwave_size = find( ttt_ripple >= window_compute_ripple_size_sharpwave_size__ms(1) & ...
            ttt_ripple <= window_compute_ripple_size_sharpwave_size__ms(2) );
        Ripple_peak_to_peak_microvolts__Individual_Ripples = max( Ripples_waveforms( :, idx_for_computing_ripples_size_and_sharpwave_size )' )' - ...
            min ( Ripples_waveforms( :, idx_for_computing_ripples_size_and_sharpwave_size )' )';
        Ripple_peak_to_peak_microvolts__mean = trimmean( Ripple_peak_to_peak_microvolts__Individual_Ripples, ...
            percent_trimmean_for_computing_mean_Ripple_size ); % Mean, computed over all individual filtered ripples : computed in a CERTAIN WINDOW AROUND t=0
        Ripple_peak_to_peak_microvolts__sem = std( Ripple_peak_to_peak_microvolts__Individual_Ripples ) ./ ... % SEM, computed over all individual filtered ripples : computed in a CERTAIN WINDOW AROUND 0
            sqrt( size(Ripples_waveforms,1) - 1 );
        Ripple_DayNumber_for_each_Individual_Ripple = ii_recording_day * ...
            ones( size( Ripple_peak_to_peak_microvolts__Individual_Ripples ) ); % Associate the current Day Number with each Individual Ripple
        Ripple_Tetrode_Depth_for_each_Individual_Ripple = Depth_of_tetrode_nominal * ...
            ones( size( Ripple_peak_to_peak_microvolts__Individual_Ripples ) ); % Associate the current Tetrode-Depth with each Individual Ripple
        % Plot:
        subplot(2,3,2); hold on ; % Plot the average Ripple Shape + Ripple Power (Abs-Hilbert)
        plot( ttt_ripple, Ripples_waveforms_mean, 'y-', 'linewidth', 0.5 ); % Ripple shape
        plot( ttt_ripple, abs( hilbert(  Ripples_waveforms_mean ) ), 'r-', 'linewidth', 0.5 ); % Ripple Power (Abs-Hilbert) 
        % Set, etc:
        title({'     Average Ripple waveform', '  Computing  Size from INDIVIDUAL ripples:', ['Peak-to-peak size (\muV) = ', num2str(Ripple_peak_to_peak_microvolts__mean,5), ' +/- ' num2str(Ripple_peak_to_peak_microvolts__sem,3)]});
        xlabel('Time (ms)');
        ylabel('Voltage (\muV)');
        set(gca, 'xtick', xticks_psth, 'xlim', xlim_psth, 'box', 'off' );
 
    
    
    % ====================================================================
    % ======= Plotting Sharpwaves (all individual sharpwaves + mean for this day): ======
    % ====================================================================
    
    % Extract Sharpwaves by looping over individual ripples, extracting RAW-LFP, and low-pass-filtering:
    
    Sharpwaves_waveforms = [];
    for ii_ripple = 1:length(Timestamps_Ripples), % Loop over Ripples
        sampling_rate_Hz = 10^6 / CSC.data.CSC_EEG_DECIMATED_SamplePeriod_microsec ; % Sampling rate of the LFP (DECIMATED)
        window_length_cut_sec = Sharpwave_extract_npoints / sampling_rate_Hz ;
        idx_data_for_Sharpwave = find( abs( Timestamps_LFP_Sleep - Timestamps_Ripples(ii_ripple) ) <= window_length_cut_sec/2 * 10^6 ); % Cut a "window" (chunk) of LFP around each Ripple
        if ( length( idx_data_for_Sharpwave ) == Sharpwave_extract_npoints ), % If the ripple/sharpwave is not too close to a session's start/end (which will truncate the "window" that I'm cutting)
            LFP_extracted =  Samples_LFP_Sleep( idx_data_for_Sharpwave );
            Wn = Sharpwave_filter_Wn_Hz / ( sampling_rate_Hz / 2); % Lowpass, Hz (Normalized Wn to half the sampling rate)
            win_sharpwave = fir1( Sharpwave_filter_N_order, Wn, 'bandpass' ); % Lowpassing, to get the Sharpwave
            %%%win_theta = win_theta / sum( win_theta ); % No need to normalize the fir1 filter by its sum!!!
            Sharpwaves_waveforms = [ Sharpwaves_waveforms ; filtfilt( win_sharpwave, 1, LFP_extracted ) ];
        end
    end % End Loop over Ripples"

        Sharpwaves_waveforms_mean = mean( Sharpwaves_waveforms );
        ttt_sharpwave = ( ( 1:length(Sharpwaves_waveforms_mean) ) - ...
            ceil( 0.5 * length(Sharpwaves_waveforms_mean) ) ) * ...
            CSC.data.CSC_EEG_DECIMATED_SamplePeriod_microsec / ...
            10^3 ; % Times (ms), with 0 centered around the ripple ; I used the Decimated srate here, since this was the srate used to extract the Sharpwave
        idx_for_computing_ripples_size_and_sharpwave_size = find( ttt_sharpwave >= window_compute_ripple_size_sharpwave_size__ms(1) & ...
            ttt_sharpwave <= window_compute_ripple_size_sharpwave_size__ms(2) );
        Sharpwave_peak_to_peak_microvolts__mean = mean( max( Sharpwaves_waveforms( :, idx_for_computing_ripples_size_and_sharpwave_size )' ) - ...
            min ( Sharpwaves_waveforms( :, idx_for_computing_ripples_size_and_sharpwave_size )' ) ); % Mean, computed over all individual filtered sharpwaves -- computed in a CERTAIN WINDOW AROUND 0
        Sharpwave_peak_to_peak_microvolts__sem = std( max( Sharpwaves_waveforms( :, idx_for_computing_ripples_size_and_sharpwave_size )' ) - ...
            min ( Sharpwaves_waveforms( :, idx_for_computing_ripples_size_and_sharpwave_size )' ) ) ./ ... % SEM, computed over all individual filtered sharpwaves -- computed in a CERTAIN WINDOW AROUND 0
            sqrt( size(Sharpwaves_waveforms,1) - 1 );
        % Plot:
        subplot(2,3,3); % Plot the average Sharpwave Shape
        plot( ttt_sharpwave, Sharpwaves_waveforms_mean, 'r-', 'linewidth', 0.5 );
        title(['Average Sharpwave waveform ;  Depth (\mum) = ', num2str(Depth_of_tetrode_nominal)]);
        xlabel('Time (ms)');
        ylabel('Voltage (\muV)');
        set(gca, 'xtick', xticks_psth, 'xlim', xlim_sharpwave, 'box', 'off' );

    
    
 % ==================================================================
    % ======== Plot the Ripple-triggered PSTH (sp/s, smoothed), and a PSTH Normalized by its Mean -- for ALL UNITS: ============
    % ==================================================================

    psth_triggered_around_ripples__all_units = zeros( length(CSC.data.filename_associated_SPIKE_files), window_size_psth_for_computing_around_ripple_ms ) ; % Initialize
    psth_triggered_around_ripples__all_units__sp2sec = zeros( length(CSC.data.filename_associated_SPIKE_files), window_size_psth_for_computing_around_ripple_ms ) ; 
    psth_triggered_around_ripples__all_units__sp2sec__smoothed = zeros( length(CSC.data.filename_associated_SPIKE_files), window_size_psth_for_computing_around_ripple_ms ) ; 
    psth_triggered_around_ripples__all_units__smoothed__norm = zeros( length(CSC.data.filename_associated_SPIKE_files), window_size_psth_for_computing_around_ripple_ms ) ; 
    Day_single_unit = zeros( 1,length(CSC.data.filename_associated_SPIKE_files)) ; % Initialize

    
    psth_triggered_around_ripples__multi_units = zeros( 1, window_size_psth_for_computing_around_ripple_ms ) ; % Initialize
    psth_triggered_around_ripples__multi_units__sp2sec = zeros( 1, window_size_psth_for_computing_around_ripple_ms ) ; 
    psth_triggered_around_ripples__multi_units__sp2sec__smoothed = zeros( 1, window_size_psth_for_computing_around_ripple_ms ) ; 
    psth_triggered_around_ripples__multi_units__smoothed__norm = zeros( 1, window_size_psth_for_computing_around_ripple_ms ) ; 

    for ii_unit = 1:length(CSC.data.filename_associated_SPIKE_files), % Loop over Units

        for ii_ripple = 1:length(Timestamps_Ripples), % Loop over Ripples
                    
            FieldSelection = [1 1 1 1 1] ; % Will read ALL the parameters, including the Samples (waveforms)
            ExtractionMode = 1 ; % Mode 1 = "Extract All"
            ExtractionModeArray = [] ; % Will read all the data
            filename_spike = CSC.data.filename_associated_SPIKE_files{ ii_unit };
        
            [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
                Nlx2MatSpike( filename_spike, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;

            window_start_time = Timestamps_Ripples(ii_ripple) - window_size_psth_for_computing_around_ripple_ms/2 * 10^3 ; % Extract spikes in a time-window around the ripple
            window_end_time = Timestamps_Ripples(ii_ripple) + window_size_psth_for_computing_around_ripple_ms/2 * 10^3 ;
            idx_spikes_extract = find( Timestamps > window_start_time & Timestamps < window_end_time );
        
            if ( ~isempty( idx_spikes_extract ) ), % If this unit HAD at all spikes in this time window, count these spikes
                Timestamps_spikes_extract = Timestamps( idx_spikes_extract );
                Timestamps_spikes_extract_shifted = Timestamps_spikes_extract - window_start_time;
                psth__1_ms_steps = ceil( Timestamps_spikes_extract_shifted / 10^3); % PSTH with 1-ms bins
                psth_triggered_around_ripples__all_units( ii_unit, psth__1_ms_steps ) = ... % Count these spikes (add 1 to the appropriate bins in the PSTH)
                    psth_triggered_around_ripples__all_units( ii_unit, psth__1_ms_steps ) + 1 ;
            end
            
     Day_single_unit(ii_unit) = ii_recording_day;

                               
        end % End loop over Ripples
        
    
        % Convert the PSTH to Sp/s and Smooth it:
    
        psth_triggered_around_ripples__all_units__sp2sec(ii_unit,:) = ...
            psth_triggered_around_ripples__all_units(ii_unit,:) * 1000 / length(Timestamps_Ripples) ; % PSTH normalized to Sp/s
     
        psth_triggered_around_ripples__all_units__sp2sec__smoothed(ii_unit,:) = ...
            smooth( psth_triggered_around_ripples__all_units__sp2sec(ii_unit,:), npoints_smooth_psth ) ; % PSTH normalized to Sp/s AND Smoothed
        
        psth_triggered_around_ripples__all_units__smoothed__norm(ii_unit,:) = ...
            psth_triggered_around_ripples__all_units__sp2sec__smoothed(ii_unit,:) / ...
            ( mean( psth_triggered_around_ripples__all_units__sp2sec__smoothed(ii_unit,:) ) + eps ); % PSTH normalized to a Mean of 1 (adding "eps" will prevent any possible divisions by 0 if mean=0)
        
    end % End loop over Units
    
     % Now do the same for the multi-unit activity
     
     for ii_ripple = 1:length(Timestamps_Ripples), % Loop over Ripples
         
         FieldSelection = [1 1 1 1 1] ; % Will read ALL the parameters, including the Samples (waveforms)
         ExtractionMode = 1 ; % Mode 1 = "Extract All"
         ExtractionModeArray = [] ; % Will read all the data
         filename_spike_multi_unit = CSC.data.file_list_in_all_threshold_crossing_un_sorted{1,1};
         
         [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
             Nlx2MatSpike( filename_spike_multi_unit, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
         
         window_start_time = Timestamps_Ripples(ii_ripple) - window_size_psth_for_computing_around_ripple_ms/2 * 10^3 ; % Extract spikes in a time-window around the ripple
         window_end_time = Timestamps_Ripples(ii_ripple) + window_size_psth_for_computing_around_ripple_ms/2 * 10^3 ;
         idx_spikes_extract = find( Timestamps > window_start_time & Timestamps < window_end_time );
         
         if ( ~isempty( idx_spikes_extract ) ), % If this unit HAD at all spikes in this time window, count these spikes
             Timestamps_spikes_extract = Timestamps( idx_spikes_extract );
             Timestamps_spikes_extract_shifted = Timestamps_spikes_extract - window_start_time;
             psth__1_ms_steps = ceil( Timestamps_spikes_extract_shifted / 10^3); % PSTH with 1-ms bins
             psth_triggered_around_ripples__multi_units( psth__1_ms_steps ) = ... % Count these spikes (add 1 to the appropriate bins in the PSTH)
                 psth_triggered_around_ripples__multi_units( psth__1_ms_steps ) + 1 ;
         end
         
     end % End loop over Ripples


        % Convert the PSTH to Sp/s and Smooth it:
    
        psth_triggered_around_ripples__multi_units__sp2sec = ...
            psth_triggered_around_ripples__multi_units * 1000 / length(Timestamps_Ripples) ; % PSTH normalized to Sp/s
     
        psth_triggered_around_ripples__multi__sp2sec__smoothed = ...
            smooth( psth_triggered_around_ripples__multi_units__sp2sec, npoints_smooth_psth ) ; % PSTH normalized to Sp/s AND Smoothed
        
        psth_triggered_around_ripples__multi__smoothed__norm = ...
            psth_triggered_around_ripples__multi_units__sp2sec__smoothed / ...
            ( mean( psth_triggered_around_ripples__multi_units__sp2sec__smoothed ) + eps ); % PSTH normalized to a Mean of 1 (adding "eps" will prevent any possible divisions by 0 if mean=0)
   
    % Plot:
    
    if ( length(CSC.data.filename_associated_SPIKE_files) >= 1 ), % If there were any single units recorded on this tetrode
        
        
%         subplot(2,3,4); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s)
            ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
                0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
%             plot( ttt_psth, psth_triggered_around_ripples__all_units__sp2sec__smoothed', 'linewidth', 1.5 );
%             title({'Ripple-triggered PSTHs, smoothed', ...
%                 ['N units = ', num2str(length(CSC.data.filename_associated_SPIKE_files))]});
%             xlabel('Time (ms)');
%             ylabel('Firing rate (Sp/s)')
%             set(gca, 'xtick', xticks_psth, 'xlim', xlim_psth, 'box', 'off' );


% % %         subplot(2,3,5); hold on ; % Plot Ripple-triggered PSTHs, smoothed (NORMALIZED by the MEAN)
% %             ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
% %                 0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
% % %             plot( ttt_psth, psth_triggered_around_ripples__all_units__smoothed__norm', 'linewidth', 1.5 );
% % %             title({'PSTHs, smoothed, normalized by the mean', ...
% % %                 ['N units = ', num2str(length(CSC.data.filename_associated_SPIKE_files))]});        
% % %             xlabel('Time (ms)');
% % %             ylabel('Firing rate (normalized)')
% % %             set(gca, 'xtick', xticks_psth, 'xlim', xlim_psth, 'box', 'off' );

%             
%         subplot(2,3,6); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s), with "ZOOM IN": to check for phase-locking
            ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
                0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
            ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
                ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
                CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
                10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
%             plot( ttt_ripple, Ripples_waveforms_mean / max(Ripples_waveforms_mean) * max(psth_triggered_around_ripples__all_units__sp2sec__smoothed(:)), ...
%                 '-', 'color', [0.5 0.5 0.5], 'linewidth', 2.5 ); % Here I Normalized the Ripples to have a peak like the maximal PSTH
%             plot( ttt_psth, psth_triggered_around_ripples__all_units__sp2sec__smoothed', 'linewidth', 1.5 );
%             title({'     Zoom In -- Ripple-triggered PSTHs, smoothed', ['Gray = Average ripple waveform  ;   ', ...
%                 'N units = ', num2str(length(CSC.data.filename_associated_SPIKE_files))]});        
%             xlabel('Time (ms)');
%             ylabel('Firing rate (Sp/s)')
%             set(gca, 'xtick', xticks_psth_ZoomIn, 'xlim', xlim_psth_ZoomIn, 'box', 'off' );

          
        subplot(2,3,5); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s), with "ZOOM IN": to check for phase-locking
            ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
                0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
            ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
                ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
                CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
                10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
            plot( ttt_ripple, Ripples_waveforms_mean / max(Ripples_waveforms_mean) * max(psth_triggered_around_ripples__all_units__sp2sec__smoothed(:)), ...
                '-', 'color', [0.5 0.5 0.5], 'linewidth', 2.5 ); % Here I Normalized the Ripples to have a peak like the maximal PSTH
            plot( ttt_psth, psth_triggered_around_ripples__all_units__sp2sec__smoothed', 'linewidth', 1.5 );
            title({'     Zoom In -- Ripple-triggered PSTHs, smoothed', ['Gray = Average ripple waveform  ;   ', ...
                'N units = ', num2str(length(CSC.data.filename_associated_SPIKE_files))]});        
            xlabel('Time (ms)');
            ylabel('Firing rate (Sp/s)')
            set(gca, 'xtick', xticks_psth_ZoomIn, 'xlim', xlim_psth_ZoomIn, 'box', 'off' );
            
            
         subplot(2,3,6); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s), with "ZOOM IN": to check for phase-locking
            ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
                0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
            ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
                ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
                CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
                10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
            plot( ttt_ripple, Ripples_waveforms_mean / max(Ripples_waveforms_mean) * max(psth_triggered_around_ripples__multi__sp2sec__smoothed), ...
                '-', 'color', [0.5 0.5 0.5], 'linewidth', 2.5 ); % Here I Normalized the Ripples to have a peak like the maximal PSTH
            plot( ttt_psth, psth_triggered_around_ripples__multi__sp2sec__smoothed', 'linewidth', 1.5 );
            title({'     Zoom In -- Ripple-triggered PSTHs, smoothed', ['Gray = Average ripple waveform  ;   ']});        
            xlabel('Time (ms)');
            ylabel('Firing rate (Sp/s)')
            set(gca, 'xtick', xticks_psth_ZoomIn, 'xlim', xlim_psth_ZoomIn, 'box', 'off' );



% 
%             
     end % End "If there were any single units recorded on this tetrode"
    
    
    
    
    % =========================================================
    % ======= Plotting Inter-Ripple-Interval Histogram: =======
    % =========================================================
        
%     subplot(4,3,6); % Plot Inter-Ripple-Interval Histogram
        InterRippleIntervals_Intervals_ms = diff( Timestamps_Ripples ) / 10^3 ; % Inter-ripples intervals (ms)
        InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms = ...
            InterRippleIntervals_Intervals_ms( find( InterRippleIntervals_Intervals_ms <= ...
            InterRipplesInterval__discard_intervals_above_this_number_ms ) ); % Inter-ripples intervals (ms), with LArge Intervals DISCARDED
        InterRippleIntervals_Hist_normalized_to_percent_occurrence = ...
            hist( InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms, InterRipplesInterval__bins_for_histogram_ms ) / ...
            length( InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms ) * ...
            100 ;
%         % Plot:
%         bar( InterRipplesInterval__bins_for_histogram_ms, InterRippleIntervals_Hist_normalized_to_percent_occurrence, ...
%             'facecolor', 'k' );
%         % Set, etc:
%         title(['Median inter-ripple interval (ms): ', num2str(median(InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms),5)]);
%         ylabel('Percent occurrence');
%         xlabel('Inter-ripple interval (ms)')
%         set( gca, 'xlim', InterRipplesInterval__xlim_ms, 'tickdir', 'out', 'box', 'off' );
%     
    
        
    % ================================================================================================
    % ======= Plotting Ripple power (Abs-Hilbert) and Computing (1)Ripple Duration @ 20-percent-height, (2) Ripple frequency: =======
    % ================================================================================================

%     subplot(2,3,5); hold on ; % Plot Ripple power (Abs-Hilbert) and Compute Ripple Duration
        ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
            ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
            CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
            10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
        idx_ripple_peak = find( ttt_ripple == 0 );
        % Plot, and compute the ripple duration:
        Ripple_Duration_ms = []; % Initialize
        Ripple_frequency = [];% Initialize
        acceptable_ripples_counter = 0; % These are only ripples which have more then 1 sinus wave
        for ii_ripple = 1:length(Timestamps_Ripples), % Loop over Ripples
            Ripple_power = abs( hilbert( Ripples_waveforms( ii_ripple, : ) ) );
%             plot( ttt_ripple, Ripple_power, 'k-', 'linewidth', 0.5 );
            idx_below_20percent_height = find( Ripple_power <= max(Ripple_power)/5 );
            idx_above_20percent_height = find( Ripple_power > max(Ripple_power)/5 );
            idx_below_20percent_height__minus__idx_ripple_peak = idx_below_20percent_height - idx_ripple_peak ;
            idx_left_crossing_of_20percent_height = max( idx_below_20percent_height__minus__idx_ripple_peak( find( idx_below_20percent_height__minus__idx_ripple_peak < 0 ) ) );
            idx_right_crossing_of_20percent_height = min( idx_below_20percent_height__minus__idx_ripple_peak( find( idx_below_20percent_height__minus__idx_ripple_peak > 0 ) ) );
            times_crossing_of_20percent_width_ms = [idx_left_crossing_of_20percent_height , idx_right_crossing_of_20percent_height] * ...
                CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / 10^3 ;
%             plot( times_crossing_of_20percent_width_ms, max(Ripple_power)/4, 'r.', 'markersize', 7 );
            Ripple_Duration_ms = [ Ripple_Duration_ms  ; diff( times_crossing_of_20percent_width_ms ) ] ;
            Current_Ripple_waveform = Ripples_waveforms( ii_ripple, : );
            Ripple_waveform_above_20percent_height = Current_Ripple_waveform(idx_above_20percent_height);
            [iHi,iLo,iCr] = findextrema(Ripple_waveform_above_20percent_height);
            if length(iHi) > 1;
                acceptable_ripples_counter = acceptable_ripples_counter + 1;
                Ripple_frequency(acceptable_ripples_counter) = (iHi(end)-iHi(1))* ...
                CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / 10^3 ;
            else end
%             figure; plot(Current_Ripple_waveform)
%             figure; plot(Ripple_waveform_above_20percent_height)
%             figure; plot(abs( hilbert( Ripples_waveforms( ii_ripple, : ) ) ))

        end % End "Loop over Ripples"
%         % Set, etc:
%         title({'Abs(Hilbert(Ripples))) + 20%-crossings (Red dots)      ', ...
%             'Ripple Duration @ 20%-height (mean +/- STD, ms):      ', ... 
%             [num2str( mean( Ripple_Duration_ms ), 3 ), ' +/- ', num2str( std( Ripple_Duration_ms ), 3 ) ] } );
%         xlabel('Time (ms)');
%         ylabel('Voltage (\muV)');
%         
    
    % ====================================================================
    % ======= Collecting some data for further Population Analysis: ======
    % ====================================================================
        
    PSD_Around_Ripples_SUM__all = [ PSD_Around_Ripples_SUM__all ; PSD_Around_Ripples_SUM ] ;
    PSD_Sleep_ALL__all = [ PSD_Sleep_ALL__all ; PSD_Sleep_ALL ] ;
    freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all = [ freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all ; freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep ];
    Ripples_waveforms_mean__all = [ Ripples_waveforms_mean__all ; Ripples_waveforms_mean ];
    psth_triggered_around_ripples__all_units__sp2sec__all = [ psth_triggered_around_ripples__all_units__sp2sec__all ; psth_triggered_around_ripples__all_units__sp2sec ];
    psth_triggered_around_ripples__all_units__sp2sec__smoothed__all = [ psth_triggered_around_ripples__all_units__sp2sec__smoothed__all ; psth_triggered_around_ripples__all_units__sp2sec__smoothed ];
    psth_triggered_around_ripples__all_units__smoothed__norm__all = [ psth_triggered_around_ripples__all_units__smoothed__norm__all ; psth_triggered_around_ripples__all_units__smoothed__norm ];
    psth_triggered_around_ripples__multi_units__sp2sec__all = [ psth_triggered_around_ripples__multi_units__sp2sec__all ; psth_triggered_around_ripples__multi_units__sp2sec ];
    psth_triggered_around_ripples__multi__sp2sec__smoothed__all = [ psth_triggered_around_ripples__multi__sp2sec__smoothed__all ; psth_triggered_around_ripples__multi__sp2sec__smoothed' ];
    psth_triggered_around_ripples__multi__smoothed__norm__all = [ psth_triggered_around_ripples__multi__smoothed__norm__all ; psth_triggered_around_ripples__multi__smoothed__norm ];
    Day_of_single_unit_all = [Day_of_single_unit_all,Day_single_unit];
    Depth_of_tetrode_nominal__all = [Depth_of_tetrode_nominal__all ; Depth_of_tetrode_nominal ] ;
    Ripple_peak_to_peak_microvolts__Individual_Ripples__all = [ Ripple_peak_to_peak_microvolts__Individual_Ripples__all ; Ripple_peak_to_peak_microvolts__Individual_Ripples ] ;
    Ripple_peak_to_peak_microvolts__mean__all = [ Ripple_peak_to_peak_microvolts__mean__all ; Ripple_peak_to_peak_microvolts__mean ] ;
    Ripple_peak_to_peak_microvolts__sem__all = [ Ripple_peak_to_peak_microvolts__sem__all ; Ripple_peak_to_peak_microvolts__sem ] ;
    Ripple_DayNumber_for_each_Individual_Ripple__all = [ Ripple_DayNumber_for_each_Individual_Ripple__all ; Ripple_DayNumber_for_each_Individual_Ripple ] ;
    Ripple_Tetrode_Depth_for_each_Individual_Ripple__all = [ Ripple_Tetrode_Depth_for_each_Individual_Ripple__all ; Ripple_Tetrode_Depth_for_each_Individual_Ripple ] ;
    Sharpwave_peak_to_peak_microvolts__mean__all = [ Sharpwave_peak_to_peak_microvolts__mean__all ; Sharpwave_peak_to_peak_microvolts__mean ] ;
    Sharpwave_peak_to_peak_microvolts__sem__all = [ Sharpwave_peak_to_peak_microvolts__sem__all ; Sharpwave_peak_to_peak_microvolts__sem ] ;
    InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms__all = [ InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms__all ; InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms] ;
    InterRippleIntervals_Hist_normalized_to_percent_occurrence__all = [ InterRippleIntervals_Hist_normalized_to_percent_occurrence__all ; InterRippleIntervals_Hist_normalized_to_percent_occurrence ] ;
    Ripple_Duration_ms__all = [ Ripple_Duration_ms__all ; Ripple_Duration_ms ] ;
    Ripple_frequency_Computed_using_local_maxima_all = [Ripple_frequency_Computed_using_local_maxima_all,Ripple_frequency];
    Timestamps_Ripples__SavedSeparatelyForEachDay{ ii_recording_day } = Timestamps_Ripples ;
    Number_of_Ripples(ii_recording_day) = N_Ripples;
    
    
    
    % =========================================
    % ======= "Print" (save) the figure: ======
    % =========================================
    
    % Find the filename by finding the backslashes ('\') in the full path name:
    idx_backslashes_in_fill_path = findstr( filename_in_list{ii_recording_day}, '\' );
    filename_fig_to_save = ...
        [ dir_save_figs, '\fig__Ripple___', ...
          filename_in_list{ii_recording_day}( idx_backslashes_in_fill_path(end)+1 : end ) ];

    % Save the figure:
    eval(['print ', filename_fig_to_save, ' -f', num2str(gcf), ' -djpeg -cmyk']);
    saveas(h,filename_fig_to_save,'fig') 

    % Close the figure, to avoid memory problems (the figure was saved already to a file):
    close gcf ;
    
    
              
end % End "Loop over recording-days"





% #################################################################
% ###########  Save some of the data into a mat-file:  ############
% #################################################################

% Defining a "Depth_microns_rescaled" variable (microns) = (Depth - mean depth) ; mean computed for each bat separately:
 
Depth_microns_rescaled__all = []; % Initialize
% Depth_microns_rescaled__all( idx_days_bat_6053 ) = Depth_of_tetrode_nominal__all( idx_days_bat_6053 ) - ...
%     mean( Depth_of_tetrode_nominal__all( idx_days_bat_6053 ) ); % Remove the mean -- Bat 6053
% Depth_microns_rescaled__all( idx_days_bat_7545 ) = Depth_of_tetrode_nominal__all( idx_days_bat_7545 ) - ...
%     mean( Depth_of_tetrode_nominal__all( idx_days_bat_7545 ) ); % Remove the mean -- Bat 7545
% Depth_microns_rescaled__all = Depth_microns_rescaled__all(:) ; % Turn into a column vector

% Population -- Data:

data_ripples.data.PSD_Around_Ripples_SUM__all = PSD_Around_Ripples_SUM__all ;
data_ripples.data.PSD_Sleep_ALL__all = PSD_Sleep_ALL__all ;
data_ripples.data.freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all = freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all ;
data_ripples.data.PSD_freqs_Ripples = PSD_freqs_Ripples ;
data_ripples.data.Ripples_waveforms_mean__all = Ripples_waveforms_mean__all ;
data_ripples.data.psth_triggered_around_ripples__all_units__sp2sec__all = psth_triggered_around_ripples__all_units__sp2sec__all ;
data_ripples.data.psth_triggered_around_ripples__all_units__sp2sec__smoothed__all = psth_triggered_around_ripples__all_units__sp2sec__smoothed__all ;
data_ripples.data.psth_triggered_around_ripples__all_units__smoothed__norm__all = psth_triggered_around_ripples__all_units__smoothed__norm__all ;
data_ripples.data.psth_triggered_around_ripples__multi_units__sp2sec__all = psth_triggered_around_ripples__multi_units__sp2sec__all ;
data_ripples.data.psth_triggered_around_ripples__multi__sp2sec__smoothed__all = psth_triggered_around_ripples__multi__sp2sec__smoothed__all ;
data_ripples.data.psth_triggered_around_ripples__multi__smoothed__norm__all = psth_triggered_around_ripples__multi__smoothed__norm__all ;
data_ripples.data.Day_of_single_unit_all = Day_of_single_unit_all;
data_ripples.data.ttt_psth_ms = ttt_psth ;
data_ripples.data.ttt_ripple_ms = ttt_ripple ;
data_ripples.data.ttt_sharpwave_ms = ttt_sharpwave ;
data_ripples.data.sampling_rate_DECIMATED = 10^6 / CSC.data.CSC_EEG_DECIMATED_SamplePeriod_microsec ;
data_ripples.data.sampling_rate_nonDECIMATED = 10^6 / CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec ;
% data_ripples.data.Depth_of_tetrode_nominal__all = Depth_of_tetrode_nominal__all ;
% data_ripples.data.Depth_microns_rescaled__all = Depth_microns_rescaled__all ;
data_ripples.data.Ripple_peak_to_peak_microvolts__Individual_Ripples__all = Ripple_peak_to_peak_microvolts__Individual_Ripples__all ;
data_ripples.data.Ripple_peak_to_peak_microvolts__mean__all = Ripple_peak_to_peak_microvolts__mean__all ;
data_ripples.data.Ripple_peak_to_peak_microvolts__sem__all = Ripple_peak_to_peak_microvolts__sem__all ;
data_ripples.data.Ripple_DayNumber_for_each_Individual_Ripple__all = Ripple_DayNumber_for_each_Individual_Ripple__all ;
data_ripples.data.Ripple_Tetrode_Depth_for_each_Individual_Ripple__all = Ripple_Tetrode_Depth_for_each_Individual_Ripple__all ;
data_ripples.data.Sharpwave_peak_to_peak_microvolts__mean__all = Sharpwave_peak_to_peak_microvolts__mean__all ;
data_ripples.data.Sharpwave_peak_to_peak_microvolts__sem__all = Sharpwave_peak_to_peak_microvolts__sem__all ; 
data_ripples.data.idx_days_bat_6053 = idx_days_bat_6053 ;
data_ripples.data.idx_days_bat_7545 = idx_days_bat_7545 ;
data_ripples.data.idx_days_first2weeks_use_for_depth_profile = idx_days_first2weeks_use_for_depth_profile ;
data_ripples.data.InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms__all = InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms__all ;
data_ripples.data.InterRippleIntervals_Hist_normalized_to_percent_occurrence__all = InterRippleIntervals_Hist_normalized_to_percent_occurrence__all ;
data_ripples.data.Ripple_Duration_ms__all = Ripple_Duration_ms__all ;
data_ripples.data.Ripple_frequency_Computed_using_local_maxima_all = Ripple_frequency_Computed_using_local_maxima_all ;
data_ripples.data.Timestamps_Ripples__SavedSeparatelyForEachDay = Timestamps_Ripples__SavedSeparatelyForEachDay ;
data_ripples.data.Number_of_Ripples = Number_of_Ripples;

% Population -- Parameters:

data_ripples.params.filename_in_list = filename_in_list ;
data_ripples.params.PSD_Ripples_npoints = PSD_Ripples_npoints ;
data_ripples.params.PSD_Ripples_freqs_50Hz_harmonics_to_interpolate_in_plots = PSD_Ripples_freqs_50Hz_harmonics_to_interpolate_in_plots ;
data_ripples.params.window_compute_ripple_size_sharpwave_size__ms = window_compute_ripple_size_sharpwave_size__ms ;
data_ripples.params.percent_trimmean_for_computing_mean_Ripple_size = percent_trimmean_for_computing_mean_Ripple_size ;
% % % data_ripples.params.bins_microns_for_computing_avg_population_ripple_size = bins_microns_for_computing_avg_population_ripple_size ;
data_ripples.params.Sharpwave_extract_npoints = Sharpwave_extract_npoints ;
data_ripples.params.Sharpwave_filter_Wn_Hz = Sharpwave_filter_Wn_Hz ;
data_ripples.params.Sharpwave_filter_N_order = Sharpwave_filter_N_order ;
data_ripples.params.window_size_psth_for_computing_around_ripple_ms = window_size_psth_for_computing_around_ripple_ms ;
data_ripples.params.npoints_smooth_psth = npoints_smooth_psth ;
data_ripples.params.InterRipplesInterval__discard_intervals_above_this_number_ms = InterRipplesInterval__discard_intervals_above_this_number_ms ;
data_ripples.params.InterRipplesInterval__bins_for_histogram_ms = InterRipplesInterval__bins_for_histogram_ms ;

% Population -- Comment:

data_ripples.COMMENT = 'Data saved by: D:\Michael\Matlab\hippo_analyze_ripples_Michael.m' ;

% Save into mat-file:

eval(['save  ', filename_mat_file_save_data, '  data_ripples']);

















% ###############################################
% ###########  Population Analysis:  ############
% ###############################################


% Plotting:

h1 = figure ; 
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[2 3 28 16]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0.2 0.2 0 0]);


subplot(2,3,1); hold on ; % Plot the PSD of Ripples -- AVERAGED ACROSS ALL RECORDING DAYS
    PSD_Around_Ripples_SUM__mean = mean( PSD_Around_Ripples_SUM__all, 1 ); % Mean along the 1'st dimension (vertical)
    PSD_Sleep_ALL__mean = mean( PSD_Sleep_ALL__all, 1 ); % Mean along the 1'st dimension (vertical)
    ylimits = log10( [ min( PSD_Sleep_ALL__mean(PSD_Ripples_xlim(2)), PSD_Around_Ripples_SUM__mean(PSD_Ripples_xlim(2)) )   max([PSD_Sleep_ALL__mean(:) ; PSD_Around_Ripples_SUM__mean(:)]) ] );
    plot( PSD_freqs_Ripples, log10( PSD_Sleep_ALL__mean ), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5 );
    plot( PSD_freqs_Ripples, log10( PSD_Around_Ripples_SUM__mean ), '-', 'color', 'y', 'linewidth', 1.5 );        
    title( { 'POPULATION ANALYSYS OF RIPPLES', ['Averages across all recordings days ,  N = ', num2str(length(filename_in_list))], ...
        'Yellow: PSD around RIPPLES, Gray:  ALL Sleep data', ...
        ['Mean +/- SEM of Ripple Freq (Hz) = peak of (Yellow-Gray) = ', num2str(mean(freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all),4), ...
        ' +/- ' num2str(std(freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all)/sqrt((length(freq_of_max_diff_between_PSD_Around_Ripples_and_PSD_sleep__all))-1),2), '     '],...
        ['Mean +/- SEM of ripple Freq (Hz) Computed using local maxima per time', num2str(mean(Ripple_frequency_Computed_using_local_maxima_all),4), ...
        ' +/- ' num2str(std(Ripple_frequency_Computed_using_local_maxima_all)/sqrt((length(Ripple_frequency_Computed_using_local_maxima_all))-1),2), '     ']} );
    xlabel( 'Frequency (Hz)' );
    ylabel('log_1_0 Power');
    set( gca, 'xlim', PSD_Ripples_xlim, 'ylim', ylimits, 'xtick', PSD_Ripples_xticks, 'box', 'off' );
 

subplot(2,3,3); hold on ; % Plot the average Ripple Shape -- AVERAGED ACROSS ALL RECORDING DAYS -- plus Ripple Power (Abs-Hilbert)
    ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
        ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
        CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
        10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
    plot( ttt_ripple, mean( Ripples_waveforms_mean__all ), 'y-', 'linewidth', 0.5 ); % Average Ripple Shape
    plot( ttt_ripple, abs( hilbert(  mean( Ripples_waveforms_mean__all ) ) ), 'r-', 'linewidth', 0.5 ); % Average Ripple Power (Abs-Hilbert) 
    title({'Population grand average', 'Average ripple waveform'});
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    set(gca, 'xtick', xticks_psth, 'xlim', xlim_psth, 'box', 'off' );



% subplot(2,3,3); hold on ; % Plot the Depth-Profile of Ripple-size and Sharpwave-size
%     ylimits = [0 450];
%     xlimits = [-200 200];
%     plot( Depth_microns_rescaled__all, ...
%         Ripple_peak_to_peak_microvolts__mean__all, 'y.', 'markersize', 9 );
%     Depth_microns_rescaled__all__first2weeks = Depth_microns_rescaled__all( idx_days_first2weeks_use_for_depth_profile ); % Data from first 2 weeks of each bat ( = stable/reliable depth data)
%     Ripple_peak_to_peak_microvolts__mean__all__first2weeks = Ripple_peak_to_peak_microvolts__mean__all( idx_days_first2weeks_use_for_depth_profile ); % Data from first 2 weeks of each bat ( = stable/reliable depth data)
%     plot( Depth_microns_rescaled__all__first2weeks, ...
%         Ripple_peak_to_peak_microvolts__mean__all__first2weeks, 'r.', 'markersize', 9 );
%     % Plot a binned-average of the rescaled-Depth:
%     Ripple_peak_to_peak_microvolts__population__by_bin = []; % Initialize
%     Depth_microns_rescaled_bin_center = [];
%     Bin_group_for_plotting_errorbars = [];
%     for ii_bin = 1:size(bins_microns_for_computing_avg_population_ripple_size,1),
%         idx_inside_bin = find( Depth_microns_rescaled__all__first2weeks >= bins_microns_for_computing_avg_population_ripple_size(ii_bin,1) & ...
%             Depth_microns_rescaled__all__first2weeks < bins_microns_for_computing_avg_population_ripple_size(ii_bin,2) ); % Data inside the Bin
%         Bin_group_for_plotting_errorbars( idx_inside_bin ) = ii_bin ; % Defining bin "groups" (for plotting errorbars)
%         Ripple_peak_to_peak_microvolts__population__by_bin__mean( ii_bin ) = ... % peak-to-peak voltage (MEAN) of data-days inside this bin
%             mean( Ripple_peak_to_peak_microvolts__mean__all__first2weeks( idx_inside_bin ) );
%         Ripple_peak_to_peak_microvolts__population__by_bin__sem( ii_bin ) = ... % peak-to-peak voltage (SEM) of data-days inside this bin
%             std( Ripple_peak_to_peak_microvolts__mean__all__first2weeks( idx_inside_bin ) ) ./ sqrt( length(idx_inside_bin) - 1 );      
%         Depth_microns_rescaled__all__bin_center( ii_bin ) = mean( bins_microns_for_computing_avg_population_ripple_size(ii_bin,:) );
%     end
%     errorbar( Depth_microns_rescaled__all__bin_center, ...
%         Ripple_peak_to_peak_microvolts__population__by_bin__mean, ...
%         Ripple_peak_to_peak_microvolts__population__by_bin__sem, ...
%         '.-', 'color', [0.5 0.5 0.5], 'linewidth', 1.0, 'markersize', 7 );
%     % one-sided t-test of (data in central bin = bin 2) versus (data in the other bins):
%     [stam, pvalue_2sided] = ttest2( Ripple_peak_to_peak_microvolts__mean__all__first2weeks( find( Bin_group_for_plotting_errorbars == 2 ) ), ...
%         Ripple_peak_to_peak_microvolts__mean__all__first2weeks( find( Bin_group_for_plotting_errorbars ~= 2 ) ), ...
%         0.05 ) ;
%     % Set, etc:
%     title({'Depth profile of Ripple size, with binned average', 'Red dots = first 2 weeks of data for each bat'...
%         '2-sided t-test of Central Bin versus other data:', ['p = ', num2str(pvalue_2sided,3)]});
%     xlabel('Depth, mormalized (\mum) = Nonimal Depth - mean depth for each bat');
%     ylabel('Peak-to-peak size (\muV)');
%     set(gca, 'xlim', xlimits, 'ylim', ylimits, 'box', 'off' );
   

subplot(2,3,4); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s) -- AVERAGED ACROSS ALL RECORDING DAYS
    ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
        0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
    grand_average_mean_firing_rate_allunits = mean( psth_triggered_around_ripples__all_units__sp2sec__smoothed__all(:) ); % Grand average of firing-rate
    plot( ttt_psth, psth_triggered_around_ripples__all_units__sp2sec__smoothed__all', 'linewidth', 1.5 );
    plot( ttt_psth, mean( psth_triggered_around_ripples__all_units__sp2sec__smoothed__all ), 'y-', 'linewidth', 2.5 );
    plot( xlim_psth, [1 1]*grand_average_mean_firing_rate_allunits, 'y:', 'linewidth', 0.5 ); % Line that denotes the Grand average of firing-rate
    title({'Population grand average', 'Ripple-triggered PSTHs, smoothed', ...
        ['N units = ', num2str(size(psth_triggered_around_ripples__all_units__sp2sec__smoothed__all,1))]});
    xlabel('Time (ms)');
    ylabel('Firing rate (Sp/s)')
    set(gca, 'xtick', xticks_psth, 'xlim', xlim_psth, 'box', 'off' );


% % % subplot(2,3,5); hold on ; % Plot Ripple-triggered PSTHs, smoothed (normalized by the MEAN) -- AVERAGED ACROSS ALL RECORDING DAYS
% % %     ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
% % %         0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
% % %     plot( ttt_psth, nanmean( psth_triggered_around_ripples__all_units__smoothed__norm__all ), 'k-', 'linewidth', 2.5 );
% % %     title({'Population grand average', 'PSTHs, smoothed, normalized by the mean', ...
% % %         ['N units = ', num2str(size(psth_triggered_around_ripples__all_units__smoothed__norm__all,1))]});        
% % %     xlabel('Time (ms)');
% % %     ylabel('Firing rate (normalized)')
% % %     set(gca, 'xtick', xticks_psth, 'xlim', xlim_psth, 'box', 'off' );

IX_days_which_had_no_ripples = find(data_ripples.data.Number_of_Ripples==0);
Valid_day_IXs = data_ripples.data.Number_of_Ripples;
Valid_day_IXs(IX_days_which_had_no_ripples) = [];
IX_ripples_more_then_threshold = find(Valid_day_IXs>minimal_num_ripples);


% Use only the days that had more then "threshold" number of ripples
InterRippleIntervals_Valid = data_ripples.data.InterRippleIntervals_Hist_normalized_to_percent_occurrence__all(IX_ripples_more_then_threshold,:);


subplot(2,3,2); hold on ; % Plot Inter-Ripple-Interval Histogram
bar( InterRipplesInterval__bins_for_histogram_ms,...
    mean( InterRippleIntervals_Valid  ), ...
    'facecolor', 'y' );
errorbar( InterRipplesInterval__bins_for_histogram_ms, ...
    mean( InterRippleIntervals_Valid  ), ...
    std( InterRippleIntervals_Valid  ) / ...
    sqrt( size(InterRippleIntervals_Valid ,1) - 1 ), ...
    'y-', 'linestyle', 'none') ;
% Set, etc:
title({['Median inter-ripple interval (ms) = ', num2str(median(InterRippleIntervals_Intervals_LargeIntervalsDiscarded_ms__all),5)], ...
    ['Ripple Duration (mean +/- SEM, ms)) = ', num2str(mean(data_ripples.data.Ripple_Duration_ms__all),3), ' +/- ', num2str(std(data_ripples.data.Ripple_Duration_ms__all),3)]});
ylabel('Percent occurrence');
xlabel('Inter-ripple interval (ms)')
set( gca, 'xlim', InterRipplesInterval__xlim_ms, 'tickdir', 'out', 'box', 'off' );



% Use only the days that had more then "threshold" number of ripples
IX_ripples_more_then_threshold_per_day = find(data_ripples.data.Number_of_Ripples>minimal_num_ripples);
IXs_single_units = [];
for ii = 1:length(IX_ripples_more_then_threshold_per_day)
     IXs_single_units = [IXs_single_units, find(Day_of_single_unit_all == IX_ripples_more_then_threshold_per_day(ii))];
end

psth_triggered_around_ripples__units__sp2sec_Valid = psth_triggered_around_ripples__all_units__sp2sec__smoothed__all(IXs_single_units,:);


subplot(2,3,5); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s), with "ZOOM IN": to check for phase-locking -- AVERAGED ACROSS ALL RECORDING DAYS
ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
    0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
    ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
    CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
    10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
grand_average_mean_firing_rate_allunits = mean( psth_triggered_around_ripples__units__sp2sec_Valid(:) ); % Grand average of firing-rate
plot( ttt_ripple, mean( Ripples_waveforms_mean__all, 1 ) / max( mean(Ripples_waveforms_mean__all,1) ) * max( mean( psth_triggered_around_ripples__units__sp2sec_Valid ) ), ...
    '-', 'color', [0.5 0.5 0.5], 'linewidth', 2.5 ); % Here I Normalized the Ripples to have a peak like the maximal PSTH
plot( xlim_psth_ZoomIn, [1 1]*grand_average_mean_firing_rate_allunits, 'y:', 'linewidth', 0.5 ); % Line that denotes the Grand average of firing-rate
plot( ttt_psth, mean( psth_triggered_around_ripples__units__sp2sec_Valid ), 'y-', 'linewidth', 2.5 ); % Plot the MEAN
plot( ttt_psth, mean( psth_triggered_around_ripples__units__sp2sec_Valid ) + ...
    std( psth_triggered_around_ripples__units__sp2sec_Valid ) / sqrt( size(psth_triggered_around_ripples__units__sp2sec_Valid,1) - 1 ) , ...
    'y-', 'linewidth', 1.0 ); % Plot the (MEAN + SEM)
plot( ttt_psth, mean( psth_triggered_around_ripples__units__sp2sec_Valid ) - ...
    std( psth_triggered_around_ripples__units__sp2sec_Valid ) / sqrt( size(psth_triggered_around_ripples__units__sp2sec_Valid,1) - 1 ) , ...
    'y-', 'linewidth', 1.0 ); % Plot the (MEAN - SEM)
legend( 'Avg Ripple', 'Avg firing rate', 'Avg PSTH', 'location', 'SouthEast' );
title({'Zoom In on PSTH + Average Ripple waveform', 'Ripple-triggered PSTHs, smoothed', ...
    ['N units = ', num2str(size(psth_triggered_around_ripples__multi__sp2sec__smoothed__all,1))]});
xlabel('Time (ms)');
ylabel('Firing rate (Sp/s)')
set(gca, 'xtick', xticks_psth_ZoomIn, 'xlim', xlim_psth_ZoomIn, 'box', 'off' );


% Use only the days that had more then "threshold" number of ripples
psth_triggered_around_ripples__multi__sp2sec_Valid = psth_triggered_around_ripples__multi__sp2sec__smoothed__all(IX_ripples_more_then_threshold,:);

subplot(2,3,6); hold on ; % Plot Ripple-triggered PSTHs, smoothed (sp/s), with "ZOOM IN": to check for phase-locking -- AVERAGED ACROSS ALL RECORDING DAYS
ttt_psth = ( 1:window_size_psth_for_computing_around_ripple_ms ) - ...
    0.5 * window_size_psth_for_computing_around_ripple_ms ; % Times (ms), with 0 centered around the ripple
ttt_ripple = ( ( 1:length(Ripples_waveforms_mean) ) - ...
    ceil( 0.5 * length(Ripples_waveforms_mean) ) ) * ...
    CSC.data.CSC_EEG_nonDECIMATEDripples_SamplePeriod_microsec / ...
    10^3 ; % Times (ms), with 0 centered around the ripple ; I used the NON-decimated srate here, since this was the srate used when cutting out the ripple waveforms
grand_average_mean_firing_rate_allunits = mean( psth_triggered_around_ripples__multi__sp2sec_Valid(:) ); % Grand average of firing-rate
plot( ttt_ripple, mean( Ripples_waveforms_mean__all, 1 ) / max( mean(Ripples_waveforms_mean__all,1) ) * max( mean( psth_triggered_around_ripples__multi__sp2sec_Valid ) ), ...
    '-', 'color', [0.5 0.5 0.5], 'linewidth', 2.5 ); % Here I Normalized the Ripples to have a peak like the maximal PSTH
plot( xlim_psth_ZoomIn, [1 1]*grand_average_mean_firing_rate_allunits, 'g:', 'linewidth', 0.5 ); % Line that denotes the Grand average of firing-rate
plot( ttt_psth, mean( psth_triggered_around_ripples__multi__sp2sec_Valid ), 'r-', 'linewidth', 2.5 ); % Plot the MEAN
plot( ttt_psth, mean( psth_triggered_around_ripples__multi__sp2sec_Valid ) + ...
    std( psth_triggered_around_ripples__multi__sp2sec_Valid ) / sqrt( size(psth_triggered_around_ripples__multi__sp2sec_Valid,1) - 1 ) , ...
    'b-', 'linewidth', 1.0 ); % Plot the (MEAN + SEM)
plot( ttt_psth, mean( psth_triggered_around_ripples__multi__sp2sec_Valid ) - ...
    std( psth_triggered_around_ripples__multi__sp2sec_Valid ) / sqrt( size(psth_triggered_around_ripples__multi__sp2sec_Valid,1) - 1 ) , ...
    'b-', 'linewidth', 1.0 ); % Plot the (MEAN - SEM)
legend( 'Avg Ripple', 'Avg firing rate', 'Avg PSTH', 'location', 'SouthEast' );
title({'Zoom In on PSTH + Average Ripple waveform', 'Ripple-triggered PSTHs, smoothed', ...
    ['N units = ', num2str(size(psth_triggered_around_ripples__multi__sp2sec_Valid,1))]});
xlabel('Time (ms)');
ylabel('Firing rate (Sp/s)')
set(gca, 'xtick', xticks_psth_ZoomIn, 'xlim', xlim_psth_ZoomIn, 'box', 'off' );


    
    
% --- "Print" (save) the Population Figure: ---

% Save the figure:

filename_population_fig_to_save = [ dir_save_figs, '\fig__Ripples___Population' ];

eval(['print ', filename_population_fig_to_save, ' -f', num2str(gcf), ' -djpeg -cmyk']);
saveas(h1,filename_population_fig_to_save,'fig') 

    





% === End: ===


cd D:\Michael\Matlab ;


% ======= THE END ===============


