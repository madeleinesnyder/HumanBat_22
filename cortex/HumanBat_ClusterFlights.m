function [flightPaths] = HumanBat_ClusterFlights(cortex_data,segmented_trajectories,cortex_ttl_seconds,dist,splines)
    
    % Uses heirarchical clustering to cluster the flights
    % Search for cortex_ttl_seconds(1)*120:end (this replaced all)

    x_mean = segmented_trajectories.trajectories_continuous(:,cortex_ttl_seconds(1)*120:end)./1000;  
    t = [1:length(cortex_data.AnalogSignals(cortex_ttl_seconds(1)*120:end,:))]./cortex_data.VideoFrameRate; % NOT SURE ABOUT THIS
    new_t=t;

    % % Flight Room references (provvisory)
    xR = +2.85; xL = -2.85; yF = 2.50;  yB = -2.50;  zT = 2.20;                 %Flight volume coordinates
    F3 = [2.56; 1.43; 0.72];    F4 = [2.56; -1.24; 0.72];                       %Feeders coordinates
    
    %% Re-extract trajectories
    % Filter and interpolate (WORKS)
    x_filt = medfilt1(x_mean,cortex_data.VideoFrameRate/2,[],2,'omitnan','truncate'); %filter after interpolating
    %x_intr = HumanBat_interpolate_nans(x_filt);%  
    x_intr = fillmissing(x_filt,'next',2,'EndValues','nearest');
    x_intr_2 = HumanBat_interpolate_nans(x_intr);
    x_spl_pre = x_intr;
    %x_spl_pre = x_intr;
    x_spl_1 = smooth(x_spl_pre(1,:),cortex_data.VideoFrameRate/2);     x_spl_2 = smooth(x_spl_pre(2,:),cortex_data.VideoFrameRate/2);      x_spl_3 = smooth(x_spl_pre(3,:),cortex_data.VideoFrameRate/2); 
    x_spl = [x_spl_1,x_spl_2,x_spl_3]';
    x_spl=x_intr;

    %% For every segment of data, look at that and the next one and see what
    % needs to be done
    all_nan_divisions = find(isnan(x_filt(1,:)));
    for i=1:length(all_nan_divisions)-1
        if all_nan_divisions(i+1)-all_nan_divisions(i) > 1
            temp_seg_strt(i) = all_nan_divisions(i)+1;
        end
    end
    [temp_seg_strt_sorted,~] = sort(unique(temp_seg_strt));

%     cc=2; 
%     for i=4:length(temp_seg_strt_sorted)-1
%         color_ = jet(7); color_idx=1; clear segment
%         segment = x_filt(:,temp_seg_strt_sorted(i):temp_seg_strt_sorted(i+cc));
%         time_between_segments = (temp_seg_strt_sorted(i+cc) - temp_seg_strt_sorted(i))/120
%         figure(); hold on; 
%         entered_nanzone = 0;
%         for j=1:length(segment)
%             if isnan(segment(1,j)) & ~entered_nanzone
%                 entered_nanzone = 1;
%                 color_idx = color_idx+1;
%             elseif ~isnan(segment(1,j)) & entered_nanzone
%                 entered_nanzone=0;
%                 scatter3(segment(1,j),segment(2,j),segment(3,j),20,color_(color_idx,:)); xlim([-2.6 2.6]); ylim([-2.6 2.6]); zlim([0 2.6]);
%             else
%                 scatter3(segment(1,j),segment(2,j),segment(3,j),20,color_(color_idx,:)); xlim([-2.6 2.6]); ylim([-2.6 2.6]); zlim([0 2.6]);
%             end
%         end
%         hold off;
%         if time_between_segments > 3
%             disp("Landing! interpolate from landing spot.");
%             %x_filt_interp = x_filt_interp(:,temp_seg_strt_sorted(i):temp_seg_strt_sorted(i+cc))
%         else 
%             disp("No landing, interpolate on a curve.")
%         end
%     end
%%


    % Temp uncomment for plotting diferent filtered data
    xx=x_intr';
    figure(); plot3(xx(:,1),xx(:,2),xx(:,3));
    figure(); plot3(xx(1:100000,1),xx(1:100000,2),xx(1:100000,3));

    %threshold based on speed
    Vx = gradient(x_spl(1,:), 1/cortex_data.VideoFrameRate);
    Vy = gradient(x_spl(2,:), 1/cortex_data.VideoFrameRate);
    Vz = gradient(x_spl(3,:), 1/cortex_data.VideoFrameRate);
    speed = sqrt(Vx.^2 + Vy.^2 + Vz.^2); %in m/s

    nonflying = find(speed < 1.5);        toofast = find(speed > 20);
    x_flying = x_spl;                     x_flying(:,[nonflying toofast]) = nan;
    batspeed = speed;                     batspeed([nonflying toofast]) = nan;
    bflying=~isnan(batspeed)';           %vector of 1s when the bat is flying
    
    % For each sample, sum up the next 1s of data(flights are longer than 1s),Code from Nick D.
    allsums = [];
    for bf = 1 : size(bflying,1)-cortex_data.VideoFrameRate
        allsums(bf) = sum(bflying(bf:bf+cortex_data.VideoFrameRate));
    end

    % Detect flight starts and stops
    [R,rLT,rUT,rLL,rUL] = risetime(allsums);    
    [F,fLT,fUT,fLL,fUL] = falltime(allsums);           
    if length(R) ~= length(F)
    fLT(length(R)) = length(allsums);
    fUT(length(R)) = length(allsums);
    F(length(R)) = F(length(F));
    end
    flight_starts = round(rLT+cortex_data.VideoFrameRate/2);
    flight_ends = round(fLT+cortex_data.VideoFrameRate/2); %... +Fs is a sistematic correction, useful
    num_flights = size(R,2);
    ref = ones(size(flight_starts));
    avg_flight_time = mean((flight_ends-flight_starts)./cortex_data.VideoFrameRate);

    %%
    % Plot 2D flight trajectories (for given amount of data)
    %data_amt = 111100; data_strt=93000;
    data_amt=length(x_mean); data_strt=1;
    plotFlightPathsAll = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    plot(x_mean(1,data_strt:data_amt),x_mean(2,data_strt:data_amt),'.');
    hold on;        rectangle('Position',[xL yB xR-xL yF-yB]);
    scatter([F3(1) F4(1)],[F3(2) F4(2)],'filled');  hold off;
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    title(['Raw flights']);
    xlabel('m'); ylabel('m');
    hold off
    
    subplot(1,2,2);
    plot(x_spl(1,data_strt:data_amt),x_spl(2,data_strt:data_amt)); hold on; plot(x_mean(1,data_strt:data_amt),x_mean(2,data_strt:data_amt),'.','MarkerSize',1);
    rectangle('Position',[xL yB xR-xL yF-yB]);  hold off;
    xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
    title(['Spline flights: ']);
    xlabel('m'); ylabel('m');
    hold off
    
    % Plot session timeline
    flight_starts = flight_starts(data_strt<flight_starts & flight_starts<data_amt);   flight_ends = flight_ends(data_strt<flight_ends & flight_ends<data_amt);  ref = ones(1,length(flight_ends));  ref_ = ones(1,length(flight_starts));
    plotFlightTimeline = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    ax1 = subplot(3,1,1);   plot(t(data_strt:data_amt),x_mean(1,data_strt:data_amt),'.');  hold on;
    plot(new_t(data_strt:data_amt),x_spl(1,data_strt:data_amt),'--');   refline(0,F3(1));    hold off;
    legend('cluster/mean','spl');   ylabel('x (m)');
    ax2 = subplot(3,1,2);   plot(new_t(data_strt:data_amt),speed(data_strt:data_amt),'.');
    hold on;    stem(new_t(flight_starts),ref_);    stem(new_t(flight_ends),ref);  hold off;
    ylabel('v (m/s)');
    ax3 = subplot(3,1,3);  % plot(t,rew_signal);
    ylabel('Rewards');
    linkaxes([ax1,ax2,ax3],'x');    xlabel('Samples');
    
    
    % Plot flights in color time order
    R = R(data_strt<flight_starts & flight_starts<data_amt);
    plotFlightPathsStartStop = figure();
    if size(R,2) > 0
        CM = jet(size(R,2));
        for nf = 1 : size(R,2)
            hold on
            plot3(x_spl(1,flight_starts(nf):flight_ends(nf)),x_spl(2,flight_starts(nf):flight_ends(nf)),x_spl(3,flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',CM(nf,:))
            hold on
            
            fstartxyz(nf,1) = x_spl(1,flight_starts(nf)); %round(nanmean(x_spl(1,flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
            fstartxyz(nf,2) = x_spl(2,flight_starts(nf)); %round(nanmean(x_spl(2,flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
            fstartxyz(nf,3) = x_spl(3,flight_starts(nf)); %round(nanmean(x_spl(3,flight_starts(nf):flight_starts(nf)+trackData.VideoFrameRate/2)));
            
            fendxyz(nf,1) = x_spl(1,flight_ends(nf)); %round(nanmean(x_spl(1,flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
            fendxyz(nf,2) = x_spl(2,flight_ends(nf)); %round(nanmean(x_spl(2,flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
            fendxyz(nf,3) = x_spl(3,flight_ends(nf)); %round(nanmean(x_spl(3,flight_ends(nf):flight_ends(nf)+trackData.VideoFrameRate/2)));
            
            scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),50,'r','filled')
            hold on
            scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),50,'k','filled')
            %pause
        end
    else
        fstartxyz(1,1) = (0);
        fstartxyz(1,2) = (0);
        fstartxyz(1,3) = (0);
        
        fendxyz(1,1) = (0);
        fendxyz(1,2) = (0);
        fendxyz(1,3) = (0);
    end
    title(['All flights start(r)/stop(b): ']);
    % modify labels for tick marks
    view(0,90)
    xlim([-3 3])
    ylim([-3 3])
    xlabel('m'); ylabel('m');
    
    hold off

    %% Clustering flights

    % Clustering params
    ds_clus = splines;                                                                %number of 3D-points/flight for clustering 
    %madeleine 25 splines, PCA+, 1m linkage
    %angelo 6 splines, PCA-, 0.7m linakge, min 5 
    pca_features = false;                                                       %if using PCA
    k_means = false;                                                            %if using k-means
    dist = dist;%2.5;                                                                 %linkage distance
    reassign = true;                                                            %re-order clusters according to density
    N_min = 3; 
    day_index=1;

    num_flights = length(flight_starts);

    %Cut out flights, downsample to ds_clus positions per flight
    all_flights = NaN(3,max(flight_ends-flight_starts),num_flights);    %3D matrix with all flights
    all_flights_ds = NaN(3,ds_clus,num_flights);                        %3D matrix with all flights(downsampled)
    
    for nf = 1 : size(all_flights,3)
        trajectory = x_spl(:,flight_starts(nf):flight_ends(nf));
        velocity = speed(:,flight_starts(nf):flight_ends(nf));
        all_flights(:,1:(flight_ends(nf)-flight_starts(nf))+1,nf) = trajectory;
        all_flights_vel(1,1:(flight_ends(nf)-flight_starts(nf)+1),nf) = velocity;
        all_flights_ds(:,:,nf) = interp1(linspace(1,3,size(trajectory,2)),trajectory',linspace(1,3,ds_clus),'spline')';
        
        %     %Uncomment if you want to see how the downsampled flights look like
            plot3(all_flights(1,:,nf),all_flights(2,:,nf),all_flights(3,:,nf),'Color','b');
            hold on;
            plot3(all_flights_ds(1,:,nf),all_flights_ds(2,:,nf),all_flights_ds(3,:,nf),'Color','r');
            hold off;
        %     w = waitforbuttonpress;
    end
    
    % Define X matrix of features for clustering (downsampled coordinates, stacked together)
    X = [all_flights_ds(1,:,:), all_flights_ds(2,:,:), all_flights_ds(3,:,:)];
    X = reshape(X,3*size(all_flights_ds,2),size(R,2));
    X = X';     %so then X = #flights x #features
    
    % If dimensionality reduction is needed
    if pca_features
        [coeff,score,latent] = pca(X);     X = score(:,1:5);
    end
    
    % k-means or hierarchical clustering (with euclidean distance and shortest linkage)
    if k_means
        n_clusters = 15;    idx = kmeans(X,n_clusters);
    else
        plotClusterDistance = figure();
        Y = pdist(X,'euclidean');   Z = linkage(Y);
        hLines = dendrogram(Z,0);  hold on;    refline(0,dist);    hold off;
        idx = cluster(Z,'Cutoff',dist,'Criterion','distance');
        title([num2str(length(unique(idx))) ' clusters: ']);
        ylim([0 10]);
    end

    % HDBSCAN (see Humanbat_HDBSCAN_Playground.m)
    %hdbscan_clusterer = HDBSCAN(X)
    
    % Create structure with flight start stop frames, id of the trajectory
    clear flight;
    flight.strt_frame = ceil(flight_starts)';
    flight.stop_frame = ceil(flight_ends)';
    flight.pos = all_flights;
    flight.vel = all_flights_vel;
    flight.id = idx;
    flight.Fs = cortex_data.VideoFrameRate;

    % Sort structure according to cluster id
    clear flight_sorted;
    [flight_sorted.id,I] = sort(flight.id);
    flight_sorted.strt_frame = flight.strt_frame(I);
    flight_sorted.stop_frame = flight.stop_frame(I);
    flight_sorted.pos = flight.pos(:,:,I);
    flight_sorted.vel = flight.vel(:,:,I);
    flight_sorted.Fs = flight.Fs;
    flight_sorted.N = size(flight_sorted.id,1);
    
    % Assign isolated clusters to cluster #flights+1
    [Ns,b] = histc(flight_sorted.id,unique(flight_sorted.id));
    flight_sorted.id(Ns(b)<N_min) = size(all_flights,3)+1;             %flight_sorted.id(Ns(b)==1) = size(all_flights,3)+1;
    id_surv_clusters = unique(flight_sorted.id);
    n_surv_clusters = size(id_surv_clusters,1);

    % Create final structure flight.clus after re-re-assignment
    clear flightPaths;
    flightPaths.id = flight_sorted.id;
    flightPaths.flight_starts_idx = flight_sorted.strt_frame';
    flightPaths.flight_ends_idx = flight_sorted.stop_frame';
    flightPaths.pos = flight_sorted.pos;
    flightPaths.vel = flight_sorted.vel;
    flightPaths.Fs = flight_sorted.Fs;
    flightPaths.N = flight_sorted.N;
    %flightPaths.day = day_index(flightPaths.flight_starts_idx);% this is the day index...
    
    for jj=1:n_surv_clusters;
        flightPaths.id(flight_sorted.id == id_surv_clusters(jj)) = jj;
    end
    id_surv_clusters = unique(flightPaths.id);
    
    %Re-assign id for convenience, if necessary
    if reassign
        new_ord = [];
        [~,new_ord] = sort(histc(flightPaths.id,id_surv_clusters(1:end-1)),'descend');
        new_ord = [new_ord; id_surv_clusters(end)];
        new_ord = circshift(new_ord,1);
        reassign_matrix =(flightPaths.id == new_ord');
        for jj=1:n_surv_clusters;
            flightPaths.id(reassign_matrix(:,jj)) = jj;
        end  
    end 
    
    %Calculate trajectory length, duration in s and interflight (take-off to take-off)
    for ii = 1:flight_sorted.N
        flightPaths.length(ii)= arclength(flightPaths.pos(1,~isnan(flightPaths.pos(1,:,ii)),ii),flightPaths.pos(2,~isnan(flightPaths.pos(2,:,ii)),ii),flightPaths.pos(3,~isnan(flightPaths.pos(3,:,ii)),ii));
        flightPaths.dur(ii) = (flightPaths.flight_ends_idx(ii)-flightPaths.flight_starts_idx(ii))./flightPaths.Fs;
    end
    %flightPaths.ifd = diff(flightPaths.flight_starts_idx)';
    
    %group the cluster ids
    for i = 1:max(flightPaths.id)
        flightPaths.clusterIndex{i} = find(flightPaths.id == i);
    end
    
    %add Tobias specific variables for future plots
    flightPaths.trajectoriesContinous = x_intr;
    flightPaths.trajectoriesSpline = x_spl;
    flightPaths.trajectoriesRaw = x_mean;
    flightPaths.batSpeed = speed';
    flightPaths.flight_starts_xyz = flightPaths.pos(:,1,:); %starting position of each flight
    flightPaths.flight_ends_xyz = zeros(length(flightPaths.pos(1,1,:)),3); %make matrix for landing position for each flight 
    for i = 1:length(flightPaths.pos(1,1,:))
        try
            flightPaths.flight_ends_xyz(i,:) = flightPaths.pos(:,find(isnan(flightPaths.pos(1,:,i)),1)-1,i); %find last xyz position
        catch
            flightPaths.flight_ends_xyz(i,:) = flightPaths.pos(:,end,i);
        end
    end
    flightPaths.flightTimeline = plotFlightTimeline;
    flightPaths.flightPathsAll = plotFlightPathsAll;
    flightPaths.flightPathsStartStop = plotFlightPathsStartStop;
    flightPaths.clusterDistance = plotClusterDistance;
    flightPaths.ds_clus = ds_clus;                                                                %number of 3D-points/flight for clustering 
    flightPaths.pca_features = pca_features;                                                       %if using PCA
    flightPaths.linkDist = dist;                                                                 %linkage distance
    flightPaths.N_min = N_min;

    %% Visualize
    
    % Force min cluster: 
    n_surv_clusters = 10;
    
    plotFlightPathsClusterEach = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    col = hsv(n_surv_clusters);
    title(['Flight clusters:']);
    for jj=1:n_surv_clusters;
        id = find(flightPaths.id==jj);
        
        subplot(3,n_surv_clusters,jj);
        avg_take_off = [];
        for ii=1:size(id,1);
            hold on;
            title(['Cluster' num2str(jj) '  (' num2str(size(id,1)) ' flights)'])
            plot3(flightPaths.pos(1,:,id(ii)),flightPaths.pos(2,:,id(ii)),flightPaths.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
            avg_take_off = [avg_take_off flightPaths.pos(:,1,id(ii))];
            hold on;
        end
        take_off = mean(avg_take_off,2);
        %textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
        
        plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
        xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
        xlabel('x');    ylabel('y');    zlabel('z');    view(2);
        hold off;
        
        subplot(3,n_surv_clusters,n_surv_clusters+jj);
        avg_take_off = [];
        for ii=1:size(id,1);
            hold on;
            title(['Cluster' num2str(jj) '  (' num2str(size(id,1)) ' flights)'])
            plot3(flightPaths.pos(1,:,id(ii)),flightPaths.pos(2,:,id(ii)),flightPaths.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
            avg_take_off = [avg_take_off flightPaths.pos(:,1,id(ii))];
            hold on;
        end
        take_off = mean(avg_take_off,2);
        %textscatter3(take_off(1),take_off(2),take_off(3),"Take-off");
        
        plot3(x_spl(1,:),x_spl(2,:),x_spl(3,:),':','Color',[0.7 0.7 0.7],'MarkerSize',0.001);
        xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
        xlabel('x');    ylabel('y');    zlabel('z');    view(0,0);
        hold off;
        
        subplot(3,n_surv_clusters,n_surv_clusters*2+jj);
        histogram(flightPaths.dur(id));
        xlim([0 15]);   xlabel('Duration(s)');  ylabel('Counts');
    end
    
    flightPaths.flightPathsClusterEach = plotFlightPathsClusterEach;

    % Temporary
    figure(); hold on; clus=2;
    id = find(flightPaths.id==clus);
    for ii=1:size(id,1);
        hold on;
        title(['Cluster' num2str(clus) '  (' num2str(size(id,1)) ' flights)']); xlim([-2.9 2.9]); ylim([-2.6 2.6]); zlim([0 2.3]);
        plot3(flightPaths.pos(1,:,id(ii)),flightPaths.pos(2,:,id(ii)),flightPaths.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
        scatter3(flightPaths.pos(1,1,id(ii)),flightPaths.pos(2,1,id(ii)),flightPaths.pos(3,1,id(ii)),50,'b','filled');
        
        avg_take_off = [avg_take_off flightPaths.pos(:,1,id(ii))];
    end

    % flightpaths.id is correct

    % Re-sort so that flight starts are in order
    clear flight_sorted_ I
    [flight_sorted_.strt_frame,I] = sort(flightPaths.flight_starts_idx);
    flight_sorted_.id = flightPaths.id(I);
    flight_sorted_.stop_frame = flightPaths.flight_ends_idx(I);
    flight_sorted_.pos = flightPaths.pos(:,:,I);
    flight_sorted_.vel = flightPaths.vel(:,:,I);
    flight_sorted_.Fs = flightPaths.Fs;
    flight_sorted_.N = size(flight_sorted_.id,1);
    flight_sorted_.length = flightPaths.length(I);
    flight_sorted_.dur = flightPaths.dur(I);
    flight_sorted_.ifd = diff(flight_sorted_.strt_frame)';
    flight_sorted_.flight_starts_xyz = flightPaths.flight_starts_xyz(:,:,I);
    flight_sorted_.flight_ends_xyz = flightPaths.flight_ends_xyz(I,:);
    
    % Create final structure flight.clus after re-re-assignment
    flightPaths.id = flight_sorted_.id;
    flightPaths.flight_starts_idx = flight_sorted_.strt_frame';
    flightPaths.flight_ends_idx = flight_sorted_.stop_frame';
    flightPaths.pos = flight_sorted_.pos;
    flightPaths.vel = flight_sorted_.vel;
    flightPaths.Fs = flight_sorted_.Fs;
    flightPaths.N = flight_sorted_.N;
    flightPaths.length = flight_sorted_.length;
    flightPaths.dur = flight_sorted_.dur; 
    flightPaths.ifd = flight_sorted_.ifd;
    flightPaths.flight_starts_xyz = flight_sorted_.flight_starts_xyz; 
    flightPaths.flight_ends_xyz = flight_sorted_.flight_ends_xyz; 

    flightPaths.id = flight_sorted_.id;
    flightPaths.flight_starts_idx_aligned = flight_sorted_.strt_frame';
    flightPaths.flight_ends_idx = flight_sorted_.stop_frame';
    flightPaths.pos = flight_sorted_.pos;
    flightPaths.vel = flight_sorted_.vel;
    flightPaths.Fs = flight_sorted_.Fs;
    flightPaths.N = flight_sorted_.N;
    flightPaths.length = flight_sorted_.length;
    flightPaths.dur = flight_sorted_.dur; 
    flightPaths.ifd = flight_sorted_.ifd;
    flightPaths.flight_starts_xyz = flight_sorted_.flight_starts_xyz; 
    flightPaths.flight_ends_xyz = flight_sorted_.flight_ends_xyz; 

    %group the cluster ids
    for i = 1:max(flightPaths.id)
        flightPaths.clusterIndex{i} = find(flightPaths.id == i);
    end

    % Temporary
    figure(); hold on; clus=3;
    id = find(flightPaths.id==clus);
    for ii=1:size(id,1);
        hold on;
        title(['Cluster' num2str(clus) '  (' num2str(size(id,1)) ' flights)']); xlim([-2.9 2.9]); ylim([-2.6 2.6]); zlim([0 2.3]);
        plot3(flightPaths.pos(1,:,id(ii)),flightPaths.pos(2,:,id(ii)),flightPaths.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col(jj,:));
        scatter3(flightPaths.pos(1,1,id(ii)),flightPaths.pos(2,1,id(ii)),flightPaths.pos(3,1,id(ii)),50,'b','filled');  
    end

    % Look at flightpaths individually if you want
    %FIQ=2;
    %figure(); plot3(flightPaths.trajectoriesContinous(1,flightPaths.flight_starts_idx(FIQ):flightPaths.flight_ends_idx(FIQ)),flightPaths.trajectoriesContinous(2,flightPaths.flight_starts_idx(FIQ):flightPaths.flight_ends_idx(FIQ)),flightPaths.trajectoriesContinous(3,flightPaths.flight_starts_idx(FIQ):flightPaths.flight_ends_idx(FIQ)))

end