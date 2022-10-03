t0 = 31*60; % sec
t1 = 31*60+10; % sec
figure;
plot(avgMarkerPos(t0*120:t1*120,1));

% Format the tracking Marker data
[Location, Location_interp] = ImBat_formatTracking(Markers);
Location_times = [1:length(AnalogSignals)];

% Segment the flights
[out] =  ImBat_SegTrajectories(Location,Location_times);

t0_ind = out.flight_starts_indx(60);
t1_ind = out.flight_ends_indx(60);
global_sample_ts_usec(t1_ind)/1e6;

figure;
tiledlayout(
scatter(avgMarkerPos(t0_ind:t1_ind,1),avgMarkerPos(t0_ind:t1_ind,2));

xpos = HumanBat_interpolate_nans(out.trajectories_continuous(1,:));
ypos = HumanBat_interpolate_nans(out.trajectories_continuous(2,:));


speed_x = diff(xpos)*120;
speed_y = diff(ypos)*120;
speed = sqrt(speed_x.^2 + speed_y.^2);
speed = medfilt1(speed,5);
top_spd = prctile(speed, 90);



figure;
scatter(xpos(t0_ind:t1_ind),ypos(t0_ind:t1_ind),5,speed(t0_ind:t1_ind), 'filled');
colorbar();
