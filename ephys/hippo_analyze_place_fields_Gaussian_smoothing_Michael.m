% Plot place fields (without spatial-view fields), and saves data into the mat-file:
% D:\Michael\Data\Expdata_Processed\Combined_data\data_place_fields.mat.

%-------------------------
% What is a place cell? 
%-------------------------
% We defined a place cell using the criterion of information per spike, with a
% value > 0.4 indicating that this cell is a place cell with a spatially
% specific firing.

%-----------------
% Michael Yartsev
%-----------------


clear all ; close all ;


% ------------- General Parameters: ----------------
num_of_cells_to_plot_on_one_figure = 3 ; % Number of cells to plot on one figure
VT_Resolution = [720 576]; % Number of pixels in VT image: [X Y]
x_bin_size_pixels = 19 ; % x-Bin size for computing place fields
y_bin_size_pixels = 19 ; % y-Bin size for computing place fields
frames_per_second = 25 ; % Will be used for normalizing Place-Field Firing-Rates to Sp/s
time_spent_minimum = 0.3 ; % Will discard pixels in which the animal spent < this minimal time (seconds)
n_point_mean_filter = 15; % The number of points we will use for the position smoothing using the mean filter. 
min_num_of_neighboring_pixels_for_inclusion_in_area_computation = 2 ; % Will use this minimal number of high-firing-rate neighboring-pixels, for including a pixel in the Area computation
threshold_factor_firing_rate_for_inclusion_in_area_computation = 0.3 ; % Will use this factor (times the maximal rate) to determine the within-place-field pixels
min_NumSpikes_for_inclusion_into_PlaceFieldSpecificity_dataset = 45 ; % I will take place-field-specificity data only from place-fields having more spikes than this number
min_NumSpikes_for_inclusion_into_PlaceFieldStability_dataset = 40 ; % I will take r_2D (stability) data only from cells having more than this number of spikes IN EACH OF THE 2 *SLEEP* SESSIONS
h = 1.5; % This is what Hafting2005 defines as the smoothing factor which is actually
%the std of the gaussian kernel.
place_cell_definition_threshold_criterion_information_per_spike = 0.5; % A value larger then this indicates that the cell is a place cell.
bin_size_for_temporal_autocorrelogram= 10; % in msec
dir_save_figs = '/home/madeleine/mnt/server2/users/KQMS/HumanBat/Figures/PlaceFields';
filename_mat_file_save_data = strcat(exp_data_path,'ephys/logger13/extracted_data/Michael_Hippo_PlaceFields.mat'); % I will save the data into this mat-file
% --------------------------------------------------

% Make Inclusion List Struct
inclusion_list.file_list_allunits{1} = '/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/220407/ephys/logger13/extracted_data/00000_20220407_TT4_SS_03.ntt';
inclusion_list.file_list_associated_VT_files_allunits = ['/home/madeleine/mnt/server2/users/KQMS/HumanBat/2992814650/processed/220407/ephys/logger13/extracted_data/SingleUnits_220516_1020/SingleUnits_220407.mat'];
inclusion_list.number_of_spikes_in_every_session_allunits = 25;
inclusion_list.is_unit_is_pyramidal_is_active_allunits = [[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1];[1 1 1]];
inclusion_list.x_bin_size_pixels_all_units = 20;
inclusion_list.y_bin_size_pixels_all_units = 20;
inclusion_list.inclusion_list_criteria_comments = [];


% ======= Extract the list of Active Neurons: =======

load D:\Michael\Data\Expdata_Processed\Combined_data\hippo_2D_data_inclusion_list_of_neurons ; % Loading Inclusion-List of Neurons
% This mat-file is where I save all the important data, including:
% * inclusion_list.file_list_allunits = Ntt file names of all my spike-sorted units.
% * inclusion_list.file_list_associated_VT_files_allunits = associated VT files (mat-files).
% * inclusion_list.number_of_spikes_in_every_session_allunits = Number of spikes in every session
% * inclusion_list.is_unit_is_pyramidal_is_active_allunits = Assignments of neurons:
%       [0 1 0] = Not an acceptable unit (e.g. too few spikes overall, or funky waveform).
%       [1 1 0] = An acceptable cell, but is NOT active during behavioral sessions.
%       [1 1 1] = An acceptable AND active cell. ONLY FOR THESE CELLS WILL I COMPUTE PLACE FIELDS AND SPATIAL-VIEW FIELDS.
% * inclusion_list.x_bin_size_pixels_all_units = x-Bin size for computing place fields
% * inclusion_list.y_bin_size_pixels_all_units = y-Bin size for computing place fields
% * inclusion_list.inclusion_list_criteria_comments = the criteria I used for assigning neurons as [1 1 0] etc.

idx_acceptable_and_active_cells = ... % List of Active Neurons
    find( inclusion_list.is_unit_is_pyramidal_is_active_allunits(:,3) == 1 );



% ======== Loop over Acceptable+Active cells: ========

for ii_cell = 1:length( idx_acceptable_and_active_cells ), % Loop over Acceptable+Active cells


    filename_spike = inclusion_list.file_list_allunits{ ...
        idx_acceptable_and_active_cells( ii_cell ) };
    
    NlxHeader_Threshold = inclusion_list.NlxHeader_Threshold( ...
        idx_acceptable_and_active_cells( ii_cell ) );

    filename_VT = inclusion_list.file_list_associated_VT_files_allunits{...
        idx_acceptable_and_active_cells( ii_cell ) };

    number_of_spikes_in_every_session = ...
        inclusion_list.number_of_spikes_in_every_session_allunits( idx_acceptable_and_active_cells( ii_cell ), : ) ;

    nominal_depth_for_this_neuron = inclusion_list.nominal_depth_all_units(idx_acceptable_and_active_cells( ii_cell )); 
    eval(['load ', filename_VT]); % Load VT (video) data

    disp( ['Analyzing cell # ', num2str( ii_cell ), ' = ', filename_spike] );


    if ( rem( ii_cell-1, num_of_cells_to_plot_on_one_figure ) == 0 ), % If I've plotted N cells already... Generate a New figure
        figure ;
        % Some WYSIWYG options:
        set(gcf,'DefaultAxesFontSize',8);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 28 26]);
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[5 0.1 0 0]);
    end


    % ==== Extract VT Informatiom (define short names for the variables): ====

    unit=1;
    % Main VT data:
    VT_Timestamps = TT_unit(unit).Timestamps ;
    VT_location_CenterOfMass = VT.location_CenterOfMass ;


    % Arena Marking:
    mean_x_NW = VT_ArenaMarking.x_NW; mean_y_NW = VT_ArenaMarking.y_NW ;
    mean_x_NE = VT_ArenaMarking.x_NE; mean_y_NE = VT_ArenaMarking.y_NE ;
    mean_x_SE = VT_ArenaMarking.x_SE; mean_y_SE = VT_ArenaMarking.y_SE ;
    mean_x_SW = VT_ArenaMarking.x_SW; mean_y_SW = VT_ArenaMarking.y_SW ;


    % ==== Read the Ntt files, and extract the data that occureed within the session: ====

    FieldSelection = [1 1 1 1 1] ; % Will read ALL the parameters, including the Samples (waveforms)
    ExtractionMode = 1 ; % Mode 1 = "Extract All"
    ExtractionModeArray = [] ; % Will read all the data

    [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
        Nlx2MatSpike( filename_spike, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;

    idx_within_session = find( Timestamps >= VT.timestamps_session_limits(1) & ...
        Timestamps <= VT.timestamps_session_limits(2) ); % Find spikes that occured within this session

    SPIKES_Timestamps = Timestamps( idx_within_session );

    % MCS
    SPIKES_Timestamps = TT_unit(unit).Timestamps;


    % ==== Smooth the tracked positions with a point mean filter : ====
    VT_location_CenterOfMass_x = smooth(VT_location_CenterOfMass(:,1),n_point_mean_filter);
    VT_location_CenterOfMass_y = smooth(VT_location_CenterOfMass(:,2),n_point_mean_filter);
    VT_location_CenterOfMass = [VT_location_CenterOfMass_x,VT_location_CenterOfMass_y];

    % ==== Compute Video Center-Of-Mass locations when Spikes occurrred (for ALL spikes): ====

    location_CenterOfMass_AtApike = zeros( length( SPIKES_Timestamps ), 2 ) + NaN ; % Initialize

    for ii_spike = 1:length( SPIKES_Timestamps ), % Loop over Spikes
        % Video frame closest to the spike:
        [stam, idx_ClosestFrame] = min( abs( SPIKES_Timestamps( ii_spike ) - VT_Timestamps ) );
        % Dealing with some edge-problems ("idx_ClosestFrame" at the edges):
        if ( idx_ClosestFrame == 1 ), idx_ClosestFrame = 2 ; end ;
        if ( idx_ClosestFrame == length(VT_Timestamps) ), idx_ClosestFrame = length(VT_Timestamps)-1 ; end ;
        % Interpolating between spike frames to find the EXACT Center-Of-Mass when Spikes occurred:
        if ( VT_Timestamps( idx_ClosestFrame ) >= SPIKES_Timestamps( ii_spike ) ),
            X_CenterOfMass_AtSpike = ...
                interp1( [ VT_Timestamps( idx_ClosestFrame - 1 )  VT_Timestamps( idx_ClosestFrame ) ], ...
                [ VT_location_CenterOfMass( idx_ClosestFrame - 1, 1 )  VT_location_CenterOfMass( idx_ClosestFrame, 1 ) ], ...
                [ SPIKES_Timestamps( ii_spike ) ], 'linear' );
            Y_CenterOfMass_AtSpike = ...
                interp1( [ VT_Timestamps( idx_ClosestFrame - 1 )  VT_Timestamps( idx_ClosestFrame ) ], ...
                [ VT_location_CenterOfMass( idx_ClosestFrame - 1, 2 )  VT_location_CenterOfMass( idx_ClosestFrame, 2 ) ], ...
                [ SPIKES_Timestamps( ii_spike ) ], 'linear' );
        elseif ( VT_Timestamps( idx_ClosestFrame ) < SPIKES_Timestamps( ii_spike ) ),
            X_CenterOfMass_AtSpike = ...
                interp1( [ VT_Timestamps( idx_ClosestFrame )  VT_Timestamps( idx_ClosestFrame + 1 ) ], ...
                [ VT_location_CenterOfMass( idx_ClosestFrame, 1 )  VT_location_CenterOfMass( idx_ClosestFrame + 1, 1 ) ], ...
                [ SPIKES_Timestamps( ii_spike ) ], 'linear' );
            Y_CenterOfMass_AtSpike = ...
                interp1( [ VT_Timestamps( idx_ClosestFrame )  VT_Timestamps( idx_ClosestFrame + 1 ) ], ...
                [ VT_location_CenterOfMass( idx_ClosestFrame, 2 )  VT_location_CenterOfMass( idx_ClosestFrame + 1, 2 ) ], ...
                [ SPIKES_Timestamps( ii_spike ) ], 'linear' );
        end
        location_CenterOfMass_AtApike( ii_spike, :) = [ X_CenterOfMass_AtSpike  Y_CenterOfMass_AtSpike] ;
    end % end of looping over spikes


    %======= Define shorter variable-names: ==========

    x_video = VT_location_CenterOfMass( :, 1 ) ;
    y_video = VT_location_CenterOfMass( :, 2 ) ;
    x_spikes = location_CenterOfMass_AtApike( :, 1 ) ;
    y_spikes = location_CenterOfMass_AtApike( :, 2 ) ;



    % ======== Compute the Raw Place Field: =======


    mat_spike_density_raw = zeros( VT_Resolution(2)/y_bin_size_pixels, VT_Resolution(1)/x_bin_size_pixels ) + NaN ; % Initialize
    mat_timespent_density_raw = zeros( VT_Resolution(2)/y_bin_size_pixels, VT_Resolution(1)/x_bin_size_pixels ) + NaN ;% Initialize
    place_field_density_raw = zeros( VT_Resolution(2)/y_bin_size_pixels, VT_Resolution(1)/x_bin_size_pixels ) + NaN ;% Initialize

    mat_spike_density_raw = spike_map;
    mat_timespent_density_raw = occMap_ts;
    place_field_density_raw = mat_spike_density_raw ./ mat_timespent_density_raw;

    ephys_binCount_norm = zeros(1,length(ephys_binCount)); occ_bin = []; occ_bin_ts=[]; 
    for i=1:length(ephys_binCount)
        %find(all(testing==[3,2],2))
        aa = find(all(coord_pairs==coord_pairs_unique(i,:),2));
        occ_bin(i) = length(aa);
        occ_bin_ts(i) = length(aa)/120;
        if occ_bin_ts(i) < 0.005
            occ_bin_ts(i)=0;
        end
        ephys_binCount_norm(i) = ephys_binCount(i)/occ_bin_ts(i);
    end

ephys_binCount_norm_sm = smoothdata(ephys_binCount_norm);
ephys_binCount_sm = smoothdata(ephys_binCount);



    for ii_x_bin = 1 : VT_Resolution(1)/x_bin_size_pixels , % Loop over x-bins
        for ii_y_bin = 1 : VT_Resolution(2)/y_bin_size_pixels , % Loop over y-bins
            % Spike Density:
            mat_spike_density_raw( ii_y_bin, ii_x_bin ) = ... % All the data
                sum( x_spikes >= 1 + x_bin_size_pixels*(ii_x_bin-1) & ...
                x_spikes <  1 + x_bin_size_pixels*(ii_x_bin) & ...
                y_spikes >= 1 + y_bin_size_pixels*(ii_y_bin-1) & ...
                y_spikes <  1 + y_bin_size_pixels*(ii_y_bin) ) ;
            % Time-Spent Density:
            mat_timespent_density_raw( ii_y_bin, ii_x_bin ) = ...
                sum( x_video >= 1 + x_bin_size_pixels*(ii_x_bin-1) & ...
                x_video <  1 + x_bin_size_pixels*(ii_x_bin) & ...
                y_video >= 1 + y_bin_size_pixels*(ii_y_bin-1) & ...
                y_video <  1 + y_bin_size_pixels*(ii_y_bin) ) ;
            % Normalize Time-Spent Density from Video-Frames-Spent to Seconds-Spent
            mat_timespent_density_raw( ii_y_bin, ii_x_bin ) = ...
                mat_timespent_density_raw( ii_y_bin, ii_x_bin ) / frames_per_second ;
            % Discard pixels in which the animal Spent less than a certain Minimal amount of time --
            % (this is computed for the "idx_include_VT" data only, usually resulting in
            % DIFFERENT pixels being discarded for the Full data):
            if ( mat_timespent_density_raw( ii_y_bin, ii_x_bin ) < time_spent_minimum ),
                mat_timespent_density_raw( ii_y_bin, ii_x_bin ) = 0 ; % Discard this time-spent-density pixel
                mat_spike_density_raw( ii_y_bin, ii_x_bin ) = 0 ; % Discard this spike-density pixel
            end
        end
    end

    % Place Field = Spike Density / Time-Spent Density :
    warning off MATLAB:divideByZero ;
    place_field_density_raw = mat_spike_density_raw ./ mat_timespent_density_raw;
    warning on MATLAB:divideByZero ;
    
    
    % ======== Compute the smoothed Place Field: =======
    
    % Create gaussian kernel: 
    sigma = h;
    hsize = 5*round(sigma)+1;
    gaussian_kernel = fspecial('gaussian',hsize,sigma);
    %         surf(gaussian_kernel); % If you want to see how your kernel looks like


    % Smoothing = convolve with gaussian kernel: 
    mat_spike_density_smoothed = imfilter(mat_spike_density_raw,gaussian_kernel);
    mat_timespent_density_smoothed = imfilter(mat_timespent_density_raw,gaussian_kernel);


    % Place Field smoothed = Spike Density smoothed / Time-Spent Density smoothed :
    warning off MATLAB:divideByZero ;
    place_field_density_smoothed = mat_spike_density_smoothed ./ mat_timespent_density_smoothed ;
    warning on MATLAB:divideByZero ;




    % ======= Compute the PF density with NaN's at unvisited location (will later be presented as white bins in the PF figure) : ==========

    % "Legalize" a bin (remove NaN) if the bat visited any of the bin's 8 closest neighbours:

    idx_timespent_density = zeros( VT_Resolution(2)/y_bin_size_pixels, VT_Resolution(1)/x_bin_size_pixels ) + NaN ; % Initialize
    for ii_x_bin = 2 : VT_Resolution(1)/x_bin_size_pixels - 1 , % Loop over x-bins, NOT INCL. CAMERA-VIEW EDGES
        for ii_y_bin = 2 : VT_Resolution(2)/y_bin_size_pixels - 1 , % Loop over y-bins, NOT INCL. CAMERA-VIEW EDGES
            matrix_3x3_of_neighbors = ...
                mat_timespent_density_raw( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
            sum_including_the_central_bin = sum(sum( matrix_3x3_of_neighbors)); % Count the matrix_3x3_of_neighbors + the central bin itself
            if ( sum_including_the_central_bin  > 0 ), % If the animal visited any of this bin's 8 neighbors (3x3 region)
                idx_timespent_density(ii_y_bin,ii_x_bin) = 1; % Put 1 in the central bin
            else  % If the animal did NOT visit any of this bin's 8 neighbors (3x3 region)
                idx_timespent_density(ii_y_bin,ii_x_bin) = 0; % Put 0 in the central bin (later we will divide by this 0 and will get NaN for the firing-rate map)
            end
        end
    end
    
   
    % Place Field = Spike Density / Time-Spent Density :
    warning off MATLAB:divideByZero ;
    place_field_density_smoothed_with_NaN = (place_field_density_smoothed.* idx_timespent_density)./idx_timespent_density;
    mat_timespent_density_smoothed_with_NaN = (mat_timespent_density_smoothed.* idx_timespent_density)./idx_timespent_density;
    warning on MATLAB:divideByZero ;
    
    idx_notNaN_PlaceField = find( ~isnan( place_field_density_smoothed_with_NaN  ) ); % Find the indexes of non-NaN bins
    idx_isNaN_PlaceField = find( isnan( place_field_density_smoothed_with_NaN  ) ); % Find the indexes of NaN bins
    
    idx_notNaN_PlaceField_un_smoothed_rate_map = find( ~isnan( place_field_density_raw ) ); % Find the indexes of non-NaN bins
    
    


    % ======= Compute the following Place-Field Parameters: ===========
    % (i) Sparsity, (ii) Coherence - currently un-used, (iii) Information_per_spike, (iv) Information_per_second, (v) Area of place-field, (vi) peak firing rate,
    % is this a place cell (based on a criterion of information per spike > 0.4) : ====


    % -----  SPARSITY  (computed for the SMOOTHED field): -----
    % Sparsity = <r_i>^2 / <r_i^2> = sum( p_i * r_i )^2 / sum( p_i * r_i^2 )
    %    Where:
    %       r_i = firing rate in bin i ;
    %       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins.
    % See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

    r_i = place_field_density_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED with NaN Place Field
    p_i = mat_timespent_density_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
        sum( mat_timespent_density_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
    r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
    sparsity = sum( p_i .* r_i )^2 / sum( p_i .* ( r_i .^2 ) ) ; % Sparsity


    % -----  COHERENCE  (computed for the UN-SMOOTHED field!!!!!!!): -----
    % Coherence = "first-order autocorrelation of the Place field =
    % Correlation (correlation coefficient) between the vector of firing rates, r_i, and the
    % average firing rates of the 8 nearest neighbors of bin i.  Formally:
    %    Coherence = corrcoef( r_i, mean(8 neighbours of r_i) )
    % See: Muller RU, Kubie JL, J. Neurosci 9, 4101-4110 (1989).

    r_i_8neighbors_avg = place_field_density_raw ;  % Initialize
    for ii_x_bin = 2 : VT_Resolution(1)/x_bin_size_pixels - 1 , % Loop over x-bins, NOT INCL. EDGES
        for ii_y_bin = 2 : VT_Resolution(2)/y_bin_size_pixels - 1 , % Loop over y-bins, NOT INCL. EDGES
            matrix_3x3_of_neighbors = ...
                place_field_density_raw( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
            num_neighbors_notNaN = sum( ~isnan( matrix_3x3_of_neighbors([1:4 6:9]) ) ); % Here I did NOT count the central bin itself
            sum_including_the_central_bin = nansum( matrix_3x3_of_neighbors(:) );
            if ( num_neighbors_notNaN > 0 ), % If there ARE non-NaN neighbors -- only then do the average
                r_i_8neighbors_avg( ii_y_bin, ii_x_bin ) = ... % Remove the central bin -- and then average
                    ( sum_including_the_central_bin - place_field_density_raw( ii_y_bin, ii_x_bin ) ) / ...
                    num_neighbors_notNaN ;
            end
        end
    end
    r_i_8neighbors_avg = r_i_8neighbors_avg( idx_notNaN_PlaceField_un_smoothed_rate_map ); % Take only bins where the Original (unsmoothed) field was not-NaN
    r_i_UNsmoothed = place_field_density_raw( idx_notNaN_PlaceField_un_smoothed_rate_map ); % Here, Use the UN-SMOOTHED Place Field
    r_i_8neighbors_avg = r_i_8neighbors_avg(:) ; r_i_UNsmoothed = r_i_UNsmoothed(:) ; % Turn these into Column vectors
    ccc = corrcoef( r_i_8neighbors_avg, r_i_UNsmoothed ); % Correlation coefficient
    rrr = ccc(2,1);
    coherence_Ztransformed = 0.5 * log( (1+rrr)/(1-rrr) ); % Fisher's z-transform of r-values


    % -----  INFORMATION PER SPIKE  (computed for the SMOOTHED field): -----
    % Information_per_spike = sum( p_i * ( r_i / r ) * log2( r_i / r ) )
    %    Where:
    %       r_i = firing rate in bin i ;
    %       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
    %       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
    % See: Skaggs WE, McNaughton BL, Wilson MA, Barnes CA, Hippocampus 6(2), 149-172 (1996).

    r_i = place_field_density_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED Place Field
    p_i = mat_timespent_density_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
        sum( mat_timespent_density_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
    r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
    % % %             r = mean( r_i ) ;
    r = sum( r_i .* p_i );
    information_per_spike = sum( p_i .* ( r_i / r ) .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)


    % -----  INFORMATION PER SECOND = INFORMATION RATE  (computed for the SMOOTHED field): -----
    % Information_per_spike = sum( p_i * r_i * log2( r_i / r ) )
    %    Where:
    %       r_i = firing rate in bin i ;
    %       p_i = occupancy of bin i = time-spent by bat in bin i / total time spent in all bins ;
    %       r = mean( r_i ) = overall mean firing rate (mean over all the pixels)
    % See: Fyhn M, Molden S, Witter MP, Moser EI, Moser M-B, Science 305, 1258-1264 (2004).

    r_i = place_field_density_smoothed_with_NaN( idx_notNaN_PlaceField ); % Use the SMOOTHED Place Field
    p_i = mat_timespent_density_smoothed_with_NaN( idx_notNaN_PlaceField ) ./ ...
        sum( mat_timespent_density_smoothed_with_NaN( idx_notNaN_PlaceField ) ) ;
    r_i = r_i(:) ; p_i = p_i(:) ; % Turn r_i and p_i into Column vectors
    % % %             r = mean( r_i ) ;
    r = sum( r_i .* p_i );
    information_per_second = sum( p_i .* r_i .* log2( ( r_i + eps ) / r ) ) ; % I added a tiny number to avoid log(0)


    % ----- AREA OF PLACE-FIELD   (computed for the SMOOTHED field): -----
    % Area = Proportion of pixels with:  Firing-Rate >= Some_Factor * max( Firing-Rate )
    %       BUT: Counting only pixelse where:
    %           At least 2/8 of the pixel's neighbors also have:  Firing-Rate >= Some_Factor * max( Firing-Rate )

    r_i = place_field_density_smoothed_with_NaN ; % Initialize ; use SMOOTHED Place Field
    idx_high_firing_rate = find( ...
        r_i >= threshold_factor_firing_rate_for_inclusion_in_area_computation * max(r_i(:)) ); % Find Pixels with high firing rate
    matrix_that_counts_num_neighbors_with_high_firing_rate = zeros( size( r_i ) ) + NaN ; % Initialize
    for ii_x_bin = 2 : VT_Resolution(1)/x_bin_size_pixels - 1 , % Loop over x-bins, NOT INCL. EDGES
        for ii_y_bin = 2 : VT_Resolution(2)/y_bin_size_pixels - 1 , % Loop over y-bins, NOT INCL. EDGES
            matrix_3x3_of_neighbors = ...
                r_i( ii_y_bin-1 : ii_y_bin+1, ii_x_bin-1 : ii_x_bin+1 ) ;
            matrix_that_counts_num_neighbors_with_high_firing_rate( ii_y_bin, ii_x_bin ) = ...
                sum( sum( matrix_3x3_of_neighbors >= ...
                threshold_factor_firing_rate_for_inclusion_in_area_computation * max(r_i(:)) ) ) ... % Counting Neighboring Pixels with high firing rate
                - 1 ; % (-1) = Subtracted the central pixel itself
        end
    end
    idx_high_firing_rate_of_neighbors = find( matrix_that_counts_num_neighbors_with_high_firing_rate >= ...
        min_num_of_neighboring_pixels_for_inclusion_in_area_computation ); % Pixels with >= N neighboring pixels with high firing rate
    idx_Pixels_WithinField = intersect( idx_high_firing_rate , ...
        idx_high_firing_rate_of_neighbors ); % Final set of "within-field" pixels = pixels with high rate *AND* >=N neighbors with high rate)
    place_field_density_smoothed_with_NaN_WithinField = zeros( size( r_i ) ); % Initialize
    place_field_density_smoothed_with_NaN_WithinField( idx_Pixels_WithinField ) = 1 ; % A "place fields" with 1's within-field and 0's elsewhere
    area = length( idx_Pixels_WithinField ) / length( idx_notNaN_PlaceField ) ;
    % Area = (Num. of pixels with high rate *AND* >=N neighbors with high rate) / (Total Num. of non-NaN pixels)


    % ----- Peak Firing Rate   (computed for the SMOOTHED field): -----

    peak_firing_rate = max(max(place_field_density_smoothed_with_NaN));
    
    % ----- Is place cell? based on a criterion of information per spike > 0.4 : -----
    is_this_a_place_cell = information_per_spike>=place_cell_definition_threshold_criterion_information_per_spike; % a vlaue of '1' indicates this IS a place cell and '0' indicates it isn't.


    % ----- Compute the un-biased autocorrelation of the smoothed rate map : -----

    % Note: % This autocorrelation can also be done using FFT via:
    % C = ifft2(fft(mat_spike_density_smoothed_Grid_cell(1:38,1:38))*fft2(fliplr(flipud(mat_spike_density_smoothed_Grid_cell(1:38,1:38)))));
    % See the follwing web-link:
    % http://www.nabble.com/2D-Autocorrelation-to25530720.html#a25530720

    %But we do it like this:
%     %First let's place zero's in the NaN location so that we can compute
%     %the autocorrelation
    place_field_density_smoothed_zero_padded = place_field_density_smoothed_with_NaN;
    place_field_density_smoothed_zero_padded(idx_isNaN_PlaceField) = 0; % We replace the NaN's with zero's.
    [NumRow,NumColumn] = size(place_field_density_smoothed_zero_padded); % Number of rows and columns in matrix A
    C = xcorr2(place_field_density_smoothed_zero_padded);
    [NumRowC,NumColumnC] = size(C); % == [NumRow+NumRow-1, NumColumn+NumColumn-1]
    basebias = zeros(NumRowC,NumColumnC);
    bri = 1:min(NumColumn);
    basebias(1,bri) = bri;
    basebias(1,NumColumnC-length(bri)+1:NumColumnC) = fliplr(bri);
    I = find(~basebias(1,:));
    if ~isempty(I),
        basebias(1,I) = max(bri)*ones(size(I));
    end %---------------------------------------------- first row complete
    for i = 2:min(NumRow); %----------------- clone & modify for rest of rows
        basebias(i,:) = i*basebias(1,:);
    end
    basebias(NumRowC-i+1:NumRowC,:) = flipud(basebias(1:i,:));
    baseline = basebias(i,:);
    I = find(~basebias(:,1));
    if ~isempty(I),
        for i = 1:length(I)
            basebias(I(i),:) = baseline;
        end
    end
    C = C./basebias;
    


    
    
    
    

    % ======= Plot graphs for Single Cells: ========


    % ----- Line-drawing of Bat locations + dot-drawing of spike locations + Lines for Arena Walls : -----

    subplot( num_of_cells_to_plot_on_one_figure, 4, ...
        rem(ii_cell-1,num_of_cells_to_plot_on_one_figure)*4 + 1 );

    hold on;

    xlimits = [0  VT_Resolution(1)];
    ylimits = [0  VT_Resolution(2)];


    % Plot Arena Walls:
    line( [mean_x_SW mean_x_NW], [mean_y_SW mean_y_NW], 'color', 'k', 'linewidth', 0.8 ); % Left wall
    line( [mean_x_SE mean_x_NE], [mean_y_SE mean_y_NE], 'color', 'k', 'linewidth', 0.8 ); % Right wall
    line( [mean_x_NW mean_x_NE], [mean_y_NW mean_y_NE], 'color', 'k', 'linewidth', 0.8 ); % Top wall
    line( [mean_x_SW mean_x_SE], [mean_y_SW mean_y_SE], 'color', 'k', 'linewidth', 0.8 ); % Bottom wall

    % Plot Bat Trajectory:
    plot( x_video, y_video, 'color', [0.5 0.5 0.5], 'linewidth', 1 ); % Bat Trajectory

    % Plot Spikes Locations:
    plot( x_spikes, y_spikes, 'r.', 'markersize', 5 ); % Spikes

    % Texts, Set, etc:
    axis ij ; axis equal ; axis tight ;
    file_name = inclusion_list.file_list_allunits(idx_acceptable_and_active_cells(ii_cell));
    if is_this_a_place_cell
    text( xlimits(1)-diff(xlimits)*0.15, ylimits(1)-diff(ylimits)*0.25, ...
        { ['Unit # ', num2str( ii_cell ),': bat # ',file_name{1:1}(45:48),', Day ',file_name{1:1}(80:81),', ',file_name{1:1}(83:85),', Cell # for today: ',file_name{1:1}(90:91), ', N spikes = ', num2str(length(x_spikes)), ...
        '  ,In 3 sessions = ', num2str( number_of_spikes_in_every_session ), ' ,binsize used [cm X cm] = ', num2str(round(x_bin_size_pixels/VT_ArenaMarking.translation_constant__pixels_per_cm)), ' ,for h value of: ', num2str(h) ], ...
        ['This is a Place Cell (information per spike > ',num2str(place_cell_definition_threshold_criterion_information_per_spike),' )' ],...
        ['(Sparsity, Coherence, Information_per_spike, Information_per_second, Area, peak_firing_rate) = '], ...
        [num2str( [sparsity, coherence_Ztransformed, information_per_spike, information_per_second, area,peak_firing_rate], 3 ) ] }, ...
        'fontname', 'helvetica', 'fontsize', 8, 'interpreter', 'none' );
    else
    text( xlimits(1)-diff(xlimits)*0.15, ylimits(1)-diff(ylimits)*0.25, ...
        { ['Unit # ', num2str( ii_cell ),': bat # ',file_name{1:1}(45:48),', Day ',file_name{1:1}(80:81),', ',file_name{1:1}(83:85),', Cell # for today: ',file_name{1:1}(90:91), ', N spikes = ', num2str(length(x_spikes)), ...
        '  ,In 3 sessions = ', num2str( number_of_spikes_in_every_session ), ' ,binsize used [cm X cm] = ', num2str(round(x_bin_size_pixels/VT_ArenaMarking.translation_constant__pixels_per_cm)), ' ,for h value of: ', num2str(h) ], ...
        ['This is NOT a Place Cell (information per spike < ',num2str(place_cell_definition_threshold_criterion_information_per_spike),' )' ],...
        ['(Sparsity, Coherence, Information_per_spike, Information_per_second, Area, peak_firing_rate) = '], ...
        [num2str( [sparsity, coherence_Ztransformed, information_per_spike, information_per_second, area,peak_firing_rate], 3 ) ] }, ...
        'fontname', 'helvetica', 'fontsize', 8, 'interpreter', 'none' );
    
    end

    set( gca, 'xlim', xlimits, 'ylim', ylimits, 'visible', 'off');



    % ----- Place Field: -----

    subplot( num_of_cells_to_plot_on_one_figure, 4, ...
        rem(ii_cell-1,num_of_cells_to_plot_on_one_figure)*4 + 2 );


    xlimits = [0  VT_Resolution(1)/x_bin_size_pixels + 1 ];
    ylimits = [0  VT_Resolution(2)/y_bin_size_pixels + 1 ];


    hhh = imagesc( place_field_density_smoothed_with_NaN );

    % Set colormap:
    colormap_map = jet ;
    colormap_map(64,:) = [1 1 1] ; % Set the highest-value pixel to White = [1 1 1] : I will use White as the Color of NaN pixels
    cdata_mat = place_field_density_smoothed_with_NaN / max(place_field_density_smoothed_with_NaN(:)) * ...
        ( size(colormap_map,1) - 1 ) ; % Scale the colors in 'cdata' to ( 1 : length(colormap) - 1 ), so that the highest value will be reserved to NaN's
    cdata_mat( idx_isNaN_PlaceField ) = size(colormap_map,1) ; % Replace NaN values with the highest values (they will be colored white)
    colormap( colormap_map );
    set(hhh, 'cdatamapping', 'direct', 'cdata', cdata_mat );

    % Labels, Set, etc:
    axis ij ; axis equal ; axis tight ;
    set( gca, 'xlim', xlimits, 'ylim', ylimits, 'visible', 'off');


    % ----- Place Field - Within-Field only (plot 1's for within-field pixels, 0's elsewhere ): -----

    subplot( num_of_cells_to_plot_on_one_figure, 4, ...
        rem(ii_cell-1,num_of_cells_to_plot_on_one_figure)*4 + 3 );

    xlimits = [0  VT_Resolution(1)/x_bin_size_pixels + 1 ];
    ylimits = [0  VT_Resolution(2)/y_bin_size_pixels + 1 ];

    hhh = imagesc( place_field_density_smoothed_with_NaN_WithinField );

    % Set colormap:
    colormap_map = jet ;
    colormap_map(64,:) = [1 1 1] ; % Set the highest-value pixel to White = [1 1 1] : I will use White as the Color of NaN pixels
    cdata_mat = place_field_density_smoothed_with_NaN_WithinField / ( max(place_field_density_smoothed_with_NaN_WithinField(:)) + eps )* ...
        ( size(colormap_map,1) - 1 ) ; % Scale the colors in 'cdata' to ( 1 : length(colormap) - 1 ), so that the highest value will be reserved to NaN's
    cdata_mat( idx_isNaN_PlaceField ) = size(colormap_map,1) ; % Replace NaN values with the highest values (they will be colored white)
    colormap( colormap_map );
    set(hhh, 'cdatamapping', 'direct', 'cdata', cdata_mat );

    % Labels, Set, etc:
    axis ij ; axis equal ; axis tight ;
    set( gca, 'xlim', xlimits, 'ylim', ylimits, 'visible', 'off' );
    
    
    % ----- Auto-Correlation Map : -----
    subplot( num_of_cells_to_plot_on_one_figure, 4, ...
        rem(ii_cell-1,num_of_cells_to_plot_on_one_figure)*4 + 4 );
    xlimits = [0  VT_Resolution(1)/x_bin_size_pixels + 1 ];
    ylimits = [0  VT_Resolution(2)/y_bin_size_pixels + 1 ];
    
    hhh = imagesc( C );
    
    % Set colormap:
    colormap_map = jet ;
    colormap_map(64,:) = [1 1 1] ; % Set the highest-value pixel to White = [1 1 1] : I will use White as the Color of NaN pixels
    cdata_mat = C / ( max(C(:)) + eps )* ...
        ( size(colormap_map,1) - 1 ) ; % Scale the colors in 'cdata' to ( 1 : length(colormap) - 1 ), so that the highest value will be reserved to NaN's
    %cdata_mat( idx_isNaN_PlaceField_Grid_cell ) = size(colormap_map,1) ; % Replace NaN values with the highest values (they will be colored white)
    colormap( colormap_map );
    set(hhh, 'cdatamapping', 'direct', 'cdata', cdata_mat );
    
    % Labels, Set, etc:
    axis ij ; axis equal ; axis tight ;
    %set( gca, 'xlim', xlimits, 'ylim', ylimits, 'visible', 'off' );






    % ========= Save the data into 2 large variables, later to be saved into a mat-file: ========

    data_place_fields{ii_cell}.COMMENT = ...
        'Data saved by the script:  D:\Michael\Matlab\hippo_analyze_place_fields_ONLY.m' ;

    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.VT_Resolution = VT_Resolution ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.x_bin_size_pixels = x_bin_size_pixels ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.y_bin_size_pixels = y_bin_size_pixels ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.frames_per_second = frames_per_second ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.time_spent_minimum = time_spent_minimum ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.min_num_of_neighboring_pixels_for_inclusion_in_area_computation = min_num_of_neighboring_pixels_for_inclusion_in_area_computation ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.threshold_factor_firing_rate_for_inclusion_in_area_computation = threshold_factor_firing_rate_for_inclusion_in_area_computation ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.min_NumSpikes_for_inclusion_into_PlaceFieldSpecificity_dataset  = min_NumSpikes_for_inclusion_into_PlaceFieldSpecificity_dataset ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.min_NumSpikes_for_inclusion_into_PlaceFieldStability_dataset = min_NumSpikes_for_inclusion_into_PlaceFieldStability_dataset ;
    data_place_fields{ii_cell}.PARAMETERS_fixed_for_all_units.bin_size_cm_used = x_bin_size_pixels*VT_ArenaMarking.translation_constant__pixels_per_cm ;

    data_place_fields{ii_cell}.Arena_Marking.x_NW = mean_x_NW; data_place_fields{ii_cell}.Arena_Marking.y_NW = mean_y_NW ;
    data_place_fields{ii_cell}.Arena_Marking.x_NE = mean_x_NE; data_place_fields{ii_cell}.Arena_Marking.y_NE = mean_y_NE ;
    data_place_fields{ii_cell}.Arena_Marking.x_SE = mean_x_SE; data_place_fields{ii_cell}.Arena_Marking.y_SE = mean_y_SE ;
    data_place_fields{ii_cell}.Arena_Marking.x_SW = mean_x_SW; data_place_fields{ii_cell}.Arena_Marking.y_SW = mean_y_SW ;


    data_place_fields{ii_cell}.filename_spike = filename_spike ;
    data_place_fields{ii_cell}.NlxHeader_Threshold = NlxHeader_Threshold;
    data_place_fields{ii_cell}.filename_VT = filename_VT ;
    data_place_fields{ii_cell}.nominal_depth_for_this_neuron = nominal_depth_for_this_neuron;
    data_place_fields{ii_cell}.N_spikes = length(x_spikes) ;
    data_place_fields{ii_cell}.SPIKES_Timestamps = SPIKES_Timestamps ;
    data_place_fields{ii_cell}.number_of_spikes_in_every_session = number_of_spikes_in_every_session ;
    data_place_fields{ii_cell}.VT_Timestamps = VT_Timestamps ;
    data_place_fields{ii_cell}.VT_location_CenterOfMass = VT_location_CenterOfMass ;
    data_place_fields{ii_cell}.x_video = x_video ;
    data_place_fields{ii_cell}.y_video = y_video ;
    data_place_fields{ii_cell}.x_spikes = x_spikes ;
    data_place_fields{ii_cell}.y_spikes = y_spikes ;
    data_place_fields{ii_cell}.mat_spike_density_raw = mat_spike_density_raw ;
    data_place_fields{ii_cell}.mat_timespent_density_raw = mat_timespent_density_raw ;
    data_place_fields{ii_cell}.mat_spike_density_smoothed = mat_spike_density_smoothed ;
    data_place_fields{ii_cell}.mat_timespent_density_smoothed = mat_timespent_density_smoothed ;
    data_place_fields{ii_cell}.place_field_density_smoothed_with_NaN = place_field_density_smoothed_with_NaN ;
    data_place_fields{ii_cell}.mat_timespent_density_smoothed_with_NaN = mat_timespent_density_smoothed_with_NaN ; % although we are currently not using this variable
    data_place_fields{ii_cell}.place_field_density_raw = place_field_density_raw ;
    data_place_fields{ii_cell}.place_field_density_smoothed = place_field_density_smoothed ;
    data_place_fields{ii_cell}.place_field_density_smoothed_with_NaN = place_field_density_smoothed_with_NaN ;
    data_place_fields{ii_cell}.PLACE_FIELD_sparsity = sparsity ;
    data_place_fields{ii_cell}.PLACE_FIELD_coherence_Ztransformed = coherence_Ztransformed ;
    data_place_fields{ii_cell}.PLACE_FIELD_information_per_spike = information_per_spike ;
    data_place_fields{ii_cell}.PLACE_FIELD_information_per_second = information_per_second ;
    data_place_fields{ii_cell}.PLACE_FIELD_area = area ;
    data_place_fields{ii_cell}.PLACE_FIELD_peak_firing_rate = peak_firing_rate ;
    data_place_fields{ii_cell}.is_this_a_place_cell = is_this_a_place_cell; 










    % ======= "Print" (save) the figure: ======

    filename_fig_to_save = ...
        [ dir_save_figs, '\fig_Unit_Plots_PlaceFieldsY_Gaussian_Smoothing', num2str(gcf),'_for_bin_size_',...
        num2str(round(x_bin_size_pixels/VT_ArenaMarking.translation_constant__pixels_per_cm)),'_for_h_value_of_', num2str(h), '.jpg' ];

    eval(['print ', filename_fig_to_save, ' -f', num2str(gcf), ' -djpeg -cmyk']);



end % Running of the cells




% ======= Save mat-file with the data: ==========

eval(['save ', num2str( filename_mat_file_save_data ), '  data_place_fields ']);


cd D:\Michael\Matlab ;











% ####################################################
% ################ Population Graphs: ################
% ####################################################
%
%
% % ----- Collect data on Measures of place-field Specificity: -----
%
PLACE_FIELD_sparsity__allunits = zeros( 1, length( data_place_fields ) ) + NaN ; % Initialize
PLACE_FIELD_coherence_Ztransformed__allunits = zeros( 1, length( data_place_fields ) ) + NaN ;
PLACE_FIELD_information_per_spike__allunits = zeros( 1, length( data_place_fields ) ) + NaN ;
PLACE_FIELD_area__allunits = zeros( 1, length( data_place_fields ) ) + NaN ;
PLACE_FIELD_Stability_r_2D = zeros( 1, length( data_place_fields ) ) + NaN ;

for ii_cell = 1 : length( data_place_fields ), % Loop over cells
    PLACE_FIELD_sparsity = [NaN NaN] ; % Initialize
    PLACE_FIELD_coherence_Ztransformed = [NaN NaN] ;
    PLACE_FIELD_information_per_spike = [NaN NaN] ;
    PLACE_FIELD_area = [NaN NaN] ;
    if( data_place_fields{ii_cell}.N_spikes >= ...
            min_NumSpikes_for_inclusion_into_PlaceFieldSpecificity_dataset ), % ... AND If the place-field had enough spikes
        PLACE_FIELD_sparsity = data_place_fields{ii_cell}.PLACE_FIELD_sparsity ;
        PLACE_FIELD_coherence_Ztransformed = data_place_fields{ii_cell}.PLACE_FIELD_coherence_Ztransformed ;
        PLACE_FIELD_information_per_spike = data_place_fields{ii_cell}.PLACE_FIELD_information_per_spike ;
        PLACE_FIELD_area = data_place_fields{ii_cell}.PLACE_FIELD_area ;
    end % End "AND If the place-field had enough spikes"
    PLACE_FIELD_sparsity__allunits( ii_cell ) = nanmean( PLACE_FIELD_sparsity ); % Mean over BOTH sessions
    PLACE_FIELD_coherence_Ztransformed__allunits( ii_cell ) = nanmean( PLACE_FIELD_coherence_Ztransformed ); % Mean over BOTH sessions
    PLACE_FIELD_information_per_spike__allunits( ii_cell ) = nanmean( PLACE_FIELD_information_per_spike ); % Mean over BOTH sessions
    PLACE_FIELD_area__allunits( ii_cell ) = nanmean( PLACE_FIELD_area ); % Mean over BOTH sessions
end % End "Loop over cells"



% % ----- Collect data on place-field Stability (r_2D between the 2 place fields): -----
%
% PLACE_FIELD_Stability_r_2D__allunits = zeros( 1, length( data_place_fields ) ) + NaN ; % Initialize
%
% for ii_cell = 1 : length( data_place_fields ), % Loop over cells
%     if ( data_place_fields{ii_cell}.number_of_spikes_in_every_session(1) >= min_NumSpikes_for_inclusion_into_PlaceFieldStability_dataset & ... % If there were >= N spikes in aLL 3 SLEEP sessions (i.e. stable recordings)
%             data_place_fields{ii_cell}.number_of_spikes_in_every_session(3) >= min_NumSpikes_for_inclusion_into_PlaceFieldStability_dataset)
%             PLACE_FIELD_Stability_r_2D__allunits( ii_cell ) = data_place_fields{ii_cell}.PLACE_FIELD_Stability_r_2D ; % Save this data
%     end % End "If there were >= N spikes in aLL 3 SLEEP sessions (i.e. stable recordings)"
% end % End "Loop over cells"



% % ----- Discard data for cells where BOTH palce-fields are Diffuse: -----
%
% idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area = []; % Initialize
% for ii_cell = 1 : length( data_place_fields ), % Loop over cells
%     if ( data_place_fields{ii_cell}.Are_both_place_fields_Diffuse__have_low_Coherence_or_0_Area == 1 ), % If Diffuse
%         idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area = ...
%             [idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area , ii_cell];
%     end % End "If Diffuse"
% end % End "Loop over cells"
%
% PLACE_FIELD_sparsity__allunits( idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area ) = []; % Discard Cell
% PLACE_FIELD_coherence_Ztransformed__allunits( idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area ) = []; % Discard Cell
% PLACE_FIELD_information_per_spike__allunits( idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area ) = []; % Discard Cell
% PLACE_FIELD_area__allunits( idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area ) = []; % Discard Cell
% PLACE_FIELD_Stability_r_2D__allunits( idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area ) = []; % Discard Cell





% ----- Plot Figure with Population Histograms: -----

figure ;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',9);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 28 10]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[5 0.1 0 0]);


n_bins = 9 ; % Number of Bins for the histograms


subplot(1,4,1);
hist( PLACE_FIELD_sparsity__allunits, n_bins );
title({['Population analysis of place-field Specificity:  N active cells = ', num2str(length( data_place_fields ))] ...
    ['Used only sessions with N spikes >= ', num2str(min_NumSpikes_for_inclusion_into_PlaceFieldSpecificity_dataset)]})
%     ['* Num of Active cells that are NOT "Place cells" = BOTH place-fields are Diffuse:  ', num2str(length(idx_both_place_fields_Are_Diffuse__have_low_Coherence_or_0_Area))]});
xlabel('Sparsity');
ylabel('# Cells');
set( gca, 'tickdir', 'out', 'box', 'off' );

subplot(1,4,2);
xxx_bins = 0.05 : 0.1 : 1.0 ;
hist( PLACE_FIELD_coherence_Ztransformed__allunits, n_bins );
xlabel('Coherence (Z-transformed)');
ylabel('# Cells');
set( gca, 'tickdir', 'out', 'box', 'off' );

subplot(1,4,3);
hist( PLACE_FIELD_area__allunits, n_bins );
xlabel('Area');
ylabel('# Cells');
set( gca, 'tickdir', 'out', 'box', 'off' );


subplot(1,4,4);
hist( PLACE_FIELD_information_per_spike__allunits, n_bins );
xlabel('Information (bits/spike)');
ylabel('# Cells');
set( gca, 'tickdir', 'out', 'box', 'off' );


% subplot(2,4,5);
% idx_Stability_not_NaN = find( ~isnan( PLACE_FIELD_Stability_r_2D__allunits ) ); % Find not-NaN values
% hist( PLACE_FIELD_Stability_r_2D__allunits( idx_Stability_not_NaN ), n_bins );
% title({'Population analysis of place-field Stability', ...
%     ['Computed only for cells in which ALL 2 SLEEP sessions had N spikes >= ', num2str(min_NumSpikes_for_inclusion_into_PlaceFieldStability_dataset)]});
% xlabel('Stability between sessions (Spatial correlation, r_2_D )');
% ylabel('# Cells');
% set( gca, 'tickdir', 'out', 'box', 'off' );



% ======= "Print" (save) the Population figure: ======

filename_fig_to_save = ...
    [ dir_save_figs, '\fig_Unit_Plots_PlaceFields_Gaussina_Smoothing_Population' ];

eval(['print ', filename_fig_to_save, ' -f', num2str(gcf), ' -djpeg -cmyk']);





% ======= THE END ===============


