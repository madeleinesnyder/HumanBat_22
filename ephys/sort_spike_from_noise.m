function  [ Accepted_spikes,  Max_r_values ] = sort_spike_from_noise( Snippets, R_threshold, Plot,Library_file)
%[ Accepted_spikes,  Max_r_values ] = sort_spike_from_noise( Input_file, R_threshold, Plot,Library_file)
% This function runs an algorithm to sort spike shapes from noise (false detection), based on
% correlations of the spike's shape with a library of acceptable spikes. For each spike, the
% correlation is computed for the tetrode channel (one out of 4) which has the largest height.
% This correlation is computed for 3 lags (-1, 0, +1), to allow for a shift of +/-1 in the
% position of the waveform's peak -- and the largest of these 3 correlations is taken, and compared
% to the varabile 'r_threshold', to determine whether to accept or reject this spike.
%
% INPUT PARAMETERS:
%    Input_file = The matfile of a tetrode signal as generated by
%    extract_logger_data.m
%    library_file = mat-file containing a matrix with the acceptable file shapes.
%    r_threshold = The algorithm rejects a waveform if it has with correlation below this threshold
%        with ALL the library's shapes.
%
% OUTPUT PARAMETERS:
%    Accepted_spikes = logical:  '1' for an accepted spike and '0' for a rejected spike.
%    Max_r_values = maximal correlation (r) value for all the spikes (maximum over all
%        the spikes in the acceptable spikes library).

if nargin<4
    PWD = pwd;
    Library_file = fullfile(PWD, 'library_of_acceptable_spike_shapes.mat');
end

if nargin<3
    Plot = 0; % set to 1  to see spikes, 0 to see none
end

if nargin<2
    R_threshold = 0.93;
end

%% Load the library of spike shapes
SS=load(Library_file, 'library_of_acceptable_spike_shapes'); % Load the library of acceptable spike shapes;

%% Load data to sort
%Tetrode_in = load(Input_file,'Snippets');
% Snippets = Tetrode_in.Snippets;
% clear Tetrode_in
NTrigger = size(Snippets,3);

%% initialize output variables
Accepted_spikes = nan(NTrigger,1);
Max_r_values = nan( NTrigger,1);

%% start to loop
% [~,FileName]=fileparts(Input_file);
% Tetrode_ID = FileName(end);
for ii_spike = 1:NTrigger % Loop over snippets that reach detection threshold
    if mod(ii_spike,1000) == 0
        %fprintf(1,'Tetrode %s Spike %d/%d -> %.2f%%\n',Tetrode_ID, ii_spike,NTrigger,ii_spike/NTrigger*100)
        fprintf(1,'Spike %d/%d -> %.2f%%\n',ii_spike,NTrigger,ii_spike/NTrigger*100)
    end
    
    % Extract the current spike shape on all 4 channels:
    Spike_shape_4channels = squeeze( Snippets(:,:,ii_spike));
    
    % Choose the channel # for which the spike has the largest height:
    [ ~, Idx_channel_max_height ] = max( max( Spike_shape_4channels ) );
    Spike_shape = Spike_shape_4channels( :, Idx_channel_max_height )' ;
    
    % Compute the correlation coefficients with all the acceptable spike shapes -- lag 0 :
    xxx_lags = 2 : 31 ; % Use CENTRAL 30 points
    ccc = corrcoef([ Spike_shape(2:end-1); SS.library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
    rrr_vec_lag_0 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
    % Compute the correlation coefficients with all the acceptable spike shapes -- lag (+1) :
    xxx_lags = 1 : 30 ; % Use FIRST 30 points (RIGHT shift of the "Acceptable Spike Shapes matrix")
    ccc = corrcoef([ Spike_shape(2:end-1) ; SS.library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
    rrr_vec_lag_plus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
    % Compute the correlation coefficients with all the acceptable spike shapes -- lag (-1) :
    xxx_lags = 3 : 32 ; % Use LAST 30 points (LEFT shift of the "Acceptable Spike Shapes matrix")
    ccc = corrcoef([ Spike_shape(2:end-1) ; SS.library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
    rrr_vec_lag_minus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
    % Save the MAXIMAL r value -- use the maximal correlation among the 3 lags (-1,0,+1):
    [Max_r_values( ii_spike ), Imax] = max(max( [ rrr_vec_lag_0;  rrr_vec_lag_plus1;  rrr_vec_lag_minus1 ] ));
    
    
    % Determine if this spike should be Accepted (set value to '1') or Rejected (set value to '0'):
    Accepted_spikes( ii_spike ) = ( Max_r_values( ii_spike )  >=  R_threshold );
    % Accept the spike shape ('1') if its correlation with ANY of the acceptable shapes
    % is >= r_threshold ; else, reject the spike ('0').
    
    if Plot
        Fig=figure();
        Fig.Position = [1400 100 1200 400];
        subplot(1,3,1)
        % plot the current spike
        if Accepted_spikes( ii_spike )
            plot(Spike_shape/max(Spike_shape), 'r-', 'LineWidth',2)
        else
            plot(Spike_shape/max(Spike_shape), 'k-', 'LineWidth',2)
        end
        text(25, 0.75, sprintf('R = %.2f',Max_r_values( ii_spike )), 'FontSize', 30);
        title('Current spike')
        subplot(1,3,2)
        % plot the library spike with the best correlation
        plot(SS.library_of_acceptable_spike_shapes(Imax,:), 'k-','LineWidth',2)
        title('Library best match')
        % plot the correlation coefficients obtained accross all spikes of
        % the library
        subplot(1,3,3)
        plot(max( [ rrr_vec_lag_0;  rrr_vec_lag_plus1;  rrr_vec_lag_minus1 ] ), 'k-', 'LineWidth',2)
        hold on
        if Accepted_spikes( ii_spike )
            text(50, 0.1, sprintf('Accepted r = %.2f',Max_r_values( ii_spike )));
            plot(Imax, Max_r_values( ii_spike ), 'r*', 'MarkerSize',20)
        else
            text(50, 0.1, sprintf('Rejected r = %.2f',Max_r_values( ii_spike )));
            plot(Imax, Max_r_values( ii_spike ), 'b*', 'MarkerSize',20)
        end
        xlabel('Library spikes')
        ylabel('correlation coefficient')
        ylim([0 1])
        title('R with all library spikes')
        hold off
        pause(1.5)
        close(Fig)
    end
end

%% save the data
% IndT = strfind(Input_file,'snippets');
% Output_file = [Input_file(1:IndT-1) 'time' Input_file(IndT+length('snippets'):end)];
% save(Output_file, 'Accepted_spikes','Max_r_values','-append');
end