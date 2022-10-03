
%% === Cluster X Human === %%
ciholas_ = 1;
Mx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToMadeleine,Flight_Group_Matrix.Clusters{clus});
Mx_flights = Flight_Group_Matrix.AllFlightsToMadeleine(Mx_flightIdx)
Kx_flightIdx = ismember(Flight_Group_Matrix.AllFlightsToKevin,Flight_Group_Matrix.Clusters{clus});
Kx_flights = Flight_Group_Matrix.AllFlightsToKevin(Kx_flightIdx)

%%
if ciholas_ == 1
%% To Madeleine vs Kevin x Cluster (clus); Takeoff 
for i=1:length(B_ephys_data.TT_unit)
    ii=i%unit;%cool_list(i);
    s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
    
    % Take all usec times of flight takeoff CHANGE THIS FOR THE
    % CATEGORIES!!!!
    to_list_m = [];
    for j=1:length(Mx_flights)
        idx = Mx_flights(j);
        to_list_m = [to_list_m,ciholas_flight_struct_resort{idx}.fstart_trimmed/120];
    end
    t_temp_takeoff_m = to_list_m'; clear to_list;
    to_list_k = [];
    for j=1:length(Kx_flights)
        idx = Kx_flights(j);
        to_list_k = [to_list_k,ciholas_flight_struct_resort{idx}.fstart_trimmed/120];
    end
    t_temp_takeoff_k = to_list_k'; clear to_list;
    cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
    cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);

    to_list_m = [];
    for j=1:length(Mx_flights)
        idx = Mx_flights(j);
        to_list_m = [to_list_m,ciholas_flight_struct_resort{idx}.fend_trimmed/120];
    end
    t_temp_landing_m = to_list_m'; clear to_list;
    to_list_k = [];
    for j=1:length(Kx_flights)
        idx = Kx_flights(j);
        to_list_k = [to_list_k,ciholas_flight_struct_resort{idx}.fend_trimmed/120];
    end
    t_temp_landing_k = to_list_k'; clear to_list;

    % Add spices
    g_id_m = ones(1,length(t_temp_takeoff_m));
    clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);
    g_id_k = ones(1,length(t_temp_takeoff_k));
    clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
    interval_ = [-3 5];

    [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v5_ciholas(s',t_temp_takeoff_m',t_temp_takeoff_k',t_temp_landing_m',t_temp_landing_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To M"),strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To K"),1,1,Mx_flights,Kx_flights,ciholas_flight_struct_resort);
end


elseif cortex_ ==1
%% To Madeleine vs Kevin x Cluster (clus); Takeoff 
for i=1:length(B_ephys_data.TT_unit)
    ii=i%unit;%cool_list(i);
    s = B_ephys_data.TT_unit(ii).AlignedTimestamps/1e6;
    
    % Take all usec times of flight takeoff CHANGE THIS FOR THE
    % CATEGORIES!!!!
    to_list_m = [];
    for j=1:length(Mx_flights)
        idx = Mx_flights(j);
        %to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
        to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fstart/120];
    end
    t_temp_takeoff_m = to_list_m'; clear to_list;
    to_list_k = [];
    for j=1:length(Kx_flights)
        idx = Kx_flights(j);
        %to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart_trimmed/120];
        to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fstart/120];
    end
    t_temp_takeoff_k = to_list_k'; clear to_list;
    cv_m = cv(Flight_Group_Matrix.AllFlightsToMadeleine);
    cv_k = cv(Flight_Group_Matrix.AllFlightsToKevin);

    to_list_m = [];
    for j=1:length(Mx_flights)
        idx = Mx_flights(j);
        to_list_m = [to_list_m,cortex_flight_struct_resort{idx}.fend/120];
    end
    t_temp_landing_m = to_list_m'; clear to_list;
    to_list_k = [];
    for j=1:length(Kx_flights)
        idx = Kx_flights(j);
        to_list_k = [to_list_k,cortex_flight_struct_resort{idx}.fend/120];
    end
    t_temp_landing_k = to_list_k'; clear to_list;

    % Add spices
    g_id_m = ones(1,length(t_temp_takeoff_m));
    clr_m = repmat(lines(1),length(t_temp_takeoff_m),1);
    g_id_k = ones(1,length(t_temp_takeoff_k));
    clr_k = repmat(lines(1),length(t_temp_takeoff_k),1);
    interval_ = [-3 5];

    [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v5_cortex(s',t_temp_takeoff_m',t_temp_takeoff_k',t_temp_landing_m',t_temp_landing_k',g_id_m,g_id_k,clr_m,clr_k,interval_,[],[],strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To M"),strcat("cortexBatUnit ",num2str(ii)," : takeoff of cortexFlights To K"),1,1,Mx_flights,Kx_flights,cortex_flight_struct_resort);
end

end