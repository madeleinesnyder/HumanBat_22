function [M_flight_struct] = HumanBat_MFlightStructure(exp_data_path,B_ephys_data,ciholas_flights,cortex_flights)

% Fields:
        % M position data
        % Ciholas bat 3d position
        % Cortex 3d position
        % Ephys (Ciholas bat)
        % Ephys (Cortex bat)
        % Audio
        % K position data
    % Measures:
        % Distance to K
        % Distanct to each bat
    
    % Load in human data
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data.mat'));

    
end

