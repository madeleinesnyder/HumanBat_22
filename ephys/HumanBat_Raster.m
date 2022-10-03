function [] = HumanBat_Raster(exp_data_path,B_ephys_data,ciholas_flights)

    %% Angelo's raster code
    % Make a raster plot for each flight all putative units
    for i=1:length(B_ephys_data.TT_unit)
        s = B_ephys_data.TT_unit(i).AlignedTimestamps/1e6;
        
        % Take all usec times of flight takeoff
        t_temp_takeoff = ciholas_flights.flights.t1;
        
        % Add spices
        g_id = ones(1,length(t_temp_takeoff));
        clr = repmat(lines(1),length(t_temp_takeoff),1);
        interval_ = [-3 5];

        [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s',t_temp_takeoff',g_id,clr,interval_,strcat("Spikes from Ciholas Bat Unit ",num2str(i)," aligned to takeoff of Ciholas flight"),1,1);
    end

end
