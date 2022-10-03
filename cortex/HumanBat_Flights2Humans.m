% A closer look at the trajectories to humans

% for easier notation
VideoFrameRate=120;
F = segmented_trajectories;
window = 100;
x_spl_1 = smoothdata(posB(1,:),'movmedian',window,'omitnan');     x_spl_2 = smoothdata(posB(2,:),'movmedian',window,'omitnan');     x_spl_3 = smoothdata(posB(3,:),'movmedian',window,'omitnan');
posB_smooth = [x_spl_1;x_spl_2;x_spl_3];

kq_ctr=0; ms_ctr=0;
for i=1:length(flights)
    if flights{i}.to_kq == 1
        kq_ctr=kq_ctr+1;
        % Plot the flight in space
        figure('name',strcat("Flight #"," ",num2str(i),". ","To KQ")); 
        subplot(1,2,1); hold on; plot3(F.trajectories_continuous(1,F.flight_starts_indx(i):F.flight_ends_indx(i)),F.trajectories_continuous(2,F.flight_starts_indx(i):F.flight_ends_indx(i)),F.trajectories_continuous(3,F.flight_starts_indx(i):F.flight_ends_indx(i)));
        plot3(posM(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'r');
        plot3(posK(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'g'); hold off;
        subplot(1,2,2); hold on;
        plot3(posB(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posB(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posB(F.flight_starts_indx(i):F.flight_ends_indx(i),3));
        plot3(posM(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'r');
        plot3(posK(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'g'); hold off;
    elseif flights{i}.to_ms == 1
        ms_ctr=ms_ctr+1;
        % Plot the flight in space
        figure('name',strcat("Flight #"," ",num2str(i),". ","To MS")); 
        subplot(1,2,1); hold on; plot3(F.trajectories_continuous(1,F.flight_starts_indx(i):F.flight_ends_indx(i)),F.trajectories_continuous(2,F.flight_starts_indx(i):F.flight_ends_indx(i)),F.trajectories_continuous(3,F.flight_starts_indx(i):F.flight_ends_indx(i)));
        plot3(posM(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'r');
        plot3(posK(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'g'); hold off;
        subplot(1,2,2); hold on;
        plot3(posB(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posB(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posB(F.flight_starts_indx(i):F.flight_ends_indx(i),3));
        plot3(posM(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'r');
        plot3(posK(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'g'); hold off;
    end
end

% Plot flight 3
i=3
figure(); hold on;
plot3(posB(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posB(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posB(F.flight_starts_indx(i):F.flight_ends_indx(i),3));
plot3(posM(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posM(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'r');
plot3(posK(F.flight_starts_indx(i):F.flight_ends_indx(i),1),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),2),posK(F.flight_starts_indx(i):F.flight_ends_indx(i),3),'g'); hold off;
    
% Plot just the z
figure(); plot(posB(F.flight_starts_indx(i):F.flight_ends_indx(i),3));


