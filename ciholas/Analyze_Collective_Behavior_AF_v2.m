%% Analyze_Collective_Behavior_AF_v2
% Script for the analysis of collective behavior (April 2022)
% Requires Processed data from Ciholas recordings and Group name (D or F)

%=== Load data
extracted_BHVfile = dir(fullfile(cd, 'Extracted_Behavior_*'));          load(extracted_BHVfile.name);
disp(['Processing ', extracted_BHVfile.name, '...']);

%=== Parameters
n_ds = 100;                                                                                                     % Downsampling factor for clustering positions
n_min = 10*Fs;                                                                                                  % Minimum number of samples to be considered a cluster (at least 10s)
n_rep = 1000;                                                                                                   % Number of repetitions for the d_th shuffling
n_pairs = length(bat_pairs);                                                                                    % Number of bat pairs
edges_dist = 10.^linspace(-3,1,100);                                                                            % Edges for distance histograms (from 1 mm to 10 m)
edges_d = {r_lim(1,1):(r_lim(1,2)-r_lim(1,1))/20:r_lim(1,2) r_lim(2,1):(r_lim(2,2)-r_lim(2,1))/20:r_lim(2,2)};  % Edges for density histogram
bat_ids = [1:n_tags]';                                                                                          % Bat identities
d_th = 0.27;                                                                                                    % Threshold distance for various classifications  
for i = 1:n_tags; for j = 1:3; custom_map(:,j,i) = linspace(1,bat_clr(i,j))'; end; end                          % Custom graded colormap
if Group_name == 'D', r_fd = [-0.5,0.05,0.45];                                                                  % Feeder position  for D bats
else, r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; end                            % Feeder positions for F bats

%=== Parameters
options.save_data = 1;                                                                                          % Save Data
options.clusterFl = 1;                                                                                          % Perform Flight Clustering
options.savefigures = 1;  fig_count = 1;                                                                        % Save Figures

%=== Create analysis folder for storing the results
if options.save_data
    analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%% Define a corrective vector for excluding flight tails with v>v_th but not classified as 'flight'

b_corr = ones(size(bflying));
for i=1:n_tags
    if f_num(i)>0
        v_diff = diff(v_abs(:,i)>v_th);
        transitions = find(v_diff);
        if v_diff(transitions(1))   ==-1,b_corr(1:transitions(1),i) = 0;    end
        if v_diff(transitions(end)) == 1,b_corr(transitions(end):end,i) = 0;end
    end
end

%% Calculate interbat-distances and percentage of time with 0,1,2,... bats flying

%=== bat_dist(:,i): distance between bat_pairs(i,:)
bat_dist = zeros(T,n_pairs);
for i = 1:n_pairs,    bat_dist(:,i) = vecnorm(r(:,:,bat_pairs(i,1))-r(:,:,bat_pairs(i,2)),2,2); end

%=== Bat distances at take-off and landing
f_nn = cell(n_tags,2);  f_bd = cell(n_tags,2);
for i=1:n_tags
    %=== All distances
    for j = 1:f_num(i)
        f_bd{i,1}(j,:)= pdist2(squeeze(r(f_smp{i,1}(j,1),:,:))',squeeze(r(f_smp{i,1}(j,1),:,i)));
        f_bd{i,2}(j,:)= pdist2(squeeze(r(f_smp{i,2}(j,1),:,:))',squeeze(r(f_smp{i,2}(j,1),:,i)));
    end
    %=== NN distances
    f_nn{i,1} = min(bat_dist(f_smp{i,1},any(bat_pairs==i,2)),[],2);
    f_nn{i,2} = min(bat_dist(f_smp{i,2},any(bat_pairs==i,2)),[],2);
end

%=== Quantify percentages of time with 0,1,2,3...bats flying
perc_time = zeros(n_tags+1,1);
for i=1:n_tags+1
    perc_time(i,:) = nnz(sum(bflying,2)== i-1)/T;
end

%% Cluster position data

%=== Pool positions from all the bats and cluster them
X = []; for i = 1:n_tags; X = [X;r(~bflying(:,i)& b_corr(:,i),:,i)]; end
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.2 0.2 0.2]);
[~,centroid_all] = Cluster3D_AF_v1(X,n_min,0.2,100);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Correct feeder position(s) 
if Group_name == 'D'
    [~,fd_idx] = min(vecnorm(r_fd-centroid_all,2,2));
    r_fd(1,:) = centroid_all(fd_idx,:);                                                  
else
    for i=1:4
    [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid_all,2,2));
    if feeder_distance<0.2,r_fd(i,:) =  centroid_all(fd_idx,:);end  % Do not correct if further than 20 cm
    end
end

%=== Single bat clusters
r_clus_id = NaN(T,n_tags);  centroid = cell(n_tags,1);  samplesIn = cell(n_tags,1); p_val_clus = zeros(n_tags,1); perc_clus = zeros(n_tags,1);  CH_idx = zeros(n_tags,19);
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0 0.33 1]);
tiledlayout(n_tags,4,'TileSpacing','tight');
for i = 1:n_tags
    X = r(~bflying(:,i)& b_corr(:,i),:,i);                                                              % Consider only stationary epochs
    [id,centroid{i,1},samplesIn{i,1},p_val_clus(i,:),CH_idx(i,:)] = Cluster3D_AF_v1(X,n_min,0.2,n_ds);  % Perform clustering
    perc_clus(i,:) = nnz(id)/length(id);                                                                % fraction of points belonging to a cluster
    r_clus_id(~bflying(:,i)& b_corr(:,i),i)=id;                                                         % Assign cluster id across entire session (nan when bat is flying)
end
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Cluster flights (f_clus = 0 if options.clusterFl =0)

alpha_clus = 1.2;         %Parameter for flight clustering
f_cls = cell(n_tags,1);     
for i = 1:n_tags
    if options.clusterFl 
        f_clus(i) = FlightClus_AF_v2(squeeze(r(:,:,i)),bflying(:,i),'Alpha',alpha_clus,'Frechet',1,'Points',10);
        f_cls{i,1} = f_clus(i).sorted_ids;
        fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
    else   
        f_clus = struct([]);
        f_cls{i,1} = 0.*f_smp{i,1};  
    end
end

%% Create FLIGHTS table features associated with flights

all_tko = t(vertcat(f_smp{:,1}));
all_lnd = t(vertcat(f_smp{:,2}));

warning('off','all');
c = 1;  FLIGHTS = table();
for i = 1:n_tags
    for j=1:f_num(i)
        FLIGHTS.id(c,:) = i;                                                        % id of the bat flying
        FLIGHTS.smp1(c,:) = f_smp{i,1}(j,1);                                        % sample takeoff
        FLIGHTS.smp2(c,:) = f_smp{i,2}(j,1);                                        % sample landing
        FLIGHTS.t1(c,:) = t(f_smp{i,1}(j,1));                                       % time takeoff
        FLIGHTS.t2(c,:) = t(f_smp{i,2}(j,1));                                       % time landing
        FLIGHTS.maxN(c,:) = max(sum(bflying(f_smp{i,1}(j,1):f_smp{i,2}(j,1),:),2)); % max number of bats flying during that flight
        FLIGHTS.r1(c,:) = r(f_smp{i,1}(j,1),:,i);                                   % position at takeoff
        FLIGHTS.r2(c,:) = r(f_smp{i,2}(j,1),:,i);                                   % position at landing
        FLIGHTS.pclus1(c,:) = r_clus_id(f_smp{i,1}(j,1)-1,i);                       % positional cluster at take off
        FLIGHTS.pclus2(c,:) = r_clus_id(f_smp{i,2}(j,1)+1,i);                       % positional cluster at landing
        FLIGHTS.fclus(c,:) = f_cls{i,1}(1,j);                                       % flight cluster
        cond_f = any(vecnorm(r(f_smp{i,2}(j,1),:,i)-r_fd,2,2)<d_th);
        cond_s = any(vecnorm(r(f_smp{i,2}(j,1),:,i)-r(f_smp{i,2}(j,1),:,bat_ids(bat_ids~=i)),2,2)<d_th);
        if      cond_f && ~cond_s
            FLIGHTS.class(c,:) = 'f';                                               % exclusive feeding flight
        elseif ~cond_f && cond_s
            FLIGHTS.class(c,:) = 's';                                               % exclusive social flight
        elseif  cond_f && cond_s
            FLIGHTS.class(c,:) = 'b';                                               % feeding + social flight
        else
            FLIGHTS.class(c,:) = 'o';                                               % other flight
        end
        
        %=== Calculate time to previous landing and time to next takeoff
        if ~isempty(min(t(f_smp{i,1}(j,1))-all_lnd(all_lnd<t(f_smp{i,1}(j,1)))))
            FLIGHTS.t_prev_lnd(c,:) = min(t(f_smp{i,1}(j,1))-all_lnd(all_lnd<t(f_smp{i,1}(j,1))));
        else
            FLIGHTS.t_prev_lnd(c,:) = Inf;
        end
        if ~isempty(min(-t(f_smp{i,2}(j,1))+all_tko(all_tko>t(f_smp{i,2}(j,1)))))
            FLIGHTS.t_next_tko(c,:) = min(-t(f_smp{i,2}(j,1))+all_tko(all_tko>t(f_smp{i,2}(j,1))));
        else
            FLIGHTS.t_next_tko(c,:) = Inf;
        end
        c = c+1;
    end
end
warning('on','all');

%% Relationships between bats: Proximity, Landing and Coupling Index

%=== Initialization of matrices
%=== Proximity and coupling indexes and associated p-values
PI = zeros(n_rep+1,n_pairs);    PI_p = zeros(n_pairs,1);
CI = zeros(n_rep+1,n_pairs);    CI_p = zeros(n_pairs,1);
%=== Landing index (forward and inverse) and associated p-values
LI_1 = zeros(n_rep+1,n_pairs);  LI_2 = zeros(n_rep+1,n_pairs); LI_p = zeros(n_pairs,2);
%=== Shuffling procedure
for i = 1:n_pairs
    %=== Observed values
    PI(1,i) = nnz(bat_dist(:,i)<d_th)/T;                                    % PI: fraction of frames < d_th
    CI(1,i) = sum(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)),1)/T;% CI: fraction of frames the bats fly together
    landing_temp = f_bd{bat_pairs(i,1),2}<d_th;                             % Temporary matrix with all the landings<d_th for bat_pairs(i,1)
    LI_1(1,i) = nnz(landing_temp(:,bat_pairs(i,2)))/f_num(bat_pairs(i,1));  % LI(1,i,1) = fraction of flights for bat_pairs(i,1) landing less than d_th from bat_pairs(i,2)
    landing_temp = f_bd{bat_pairs(i,2),2}<d_th;                             % Temporary matrix with all the landings<d_th for bat_pairs(i,2)
    LI_2(1,i) = nnz(landing_temp(:,bat_pairs(i,1)))/f_num(bat_pairs(i,2));  % LI(1,i,2) = fraction of flights for bat_pairs(i,2) landing less than d_th from bat_pairs(i,1)
    frames_to_shift = randi([60*Fs T-60*Fs],1,n_rep);                       % Random number of frames between 1 min and session time-1min
    %=== Values after shuffling position and velocity
    parfor n = 1:n_rep
        r_shfl = circshift(r(:,:,bat_pairs(i,2)),frames_to_shift(n),1);     % Shift the position of the second bat in the couple
        PI(n+1,i) = nnz(vecnorm(r(:,:,bat_pairs(i,1))-r_shfl,2,2)<d_th)/T;  % Recalculate the PI
        LI_1(n+1,i) = nnz(vecnorm(r(f_smp{bat_pairs(i,1),2},:,bat_pairs(i,1))-r_shfl(f_smp{bat_pairs(i,1),2},:),2,2)<d_th)/f_num(bat_pairs(i,1)); % Calculate LI as the fraction of flights for bat_pairs(i,1) landing less than d_th from r_shfl
        r_shfl = circshift(r(:,:,bat_pairs(i,1)),frames_to_shift(n),1);     % Shift the position of the first bat in the couple
        LI_2(n+1,i) = nnz(vecnorm(r(f_smp{bat_pairs(i,2),2},:,bat_pairs(i,2))-r_shfl(f_smp{bat_pairs(i,2),2},:),2,2)<d_th)/f_num(bat_pairs(i,2)); % Calculate LI as the fraction of flights for bat_pairs(i,2) landing less than d_th from r_shfl
        b_shfl = circshift(bflying(:,bat_pairs(i,2)),frames_to_shift(n),1); % Shift the flight vector of the second bat in the couple
        CI(n+1,i) = sum(bflying(:,bat_pairs(i,1)).*b_shfl,1)/T;             % Recalculate the CI
    end
    PI_p(i,:) = nnz(PI(2:end,i)>PI(1,i))/n_rep;                                 % Calculate the p value as the fraction of PI > observed PI
    CI_p(i,:) = nnz(CI(2:end,i)>CI(1,i))/n_rep;                                 % Calculate the p value as the fraction of CI > observed CI
    LI_p(i,1) = nnz(LI_1(2:end,i)>LI_1(1,i))/n_rep;                             % Calculate the p value as the fraction of LI > observed LI
    LI_p(i,2) = nnz(LI_2(2:end,i)>LI_2(1,i))/n_rep;                             % Calculate the p value as the fraction of LI > observed LI
end

LI = cat(3,LI_1,LI_2);

%% Relationships between bats: Preferred Locations

%=== Extract first 3 preferred locations for each bat
pref_location = NaN(n_tags,3,3);    pref_time_loc = NaN(n_tags,1,3);
for i=1:n_tags
    [~,indx] = sort(samplesIn{i,1},'descend','MissingPlacement','last');                     % Sort centroids occupancy
    num_pref = min(3,numel(indx));                                                           % Check how many centroids survive
    pref_location(i,:,1:num_pref) = centroid{i,1}(indx(1:num_pref),:)';                      % Assign the first 3 preferred locations
    pref_time_loc(i,:,1:num_pref) = samplesIn{i,1}(indx(1:num_pref),:)'./T;                  % Fraction of time in the preferred locations
end

%% Relationships between bats: Coordinated flight

%=== Probability of flight (one bat vs two bats)
p_fly1 = sum(bflying,1)./T; p_fly2 = zeros(n_pairs,2);
for i = 1:length(bat_pairs)
    p_fly2(i,1) = sum(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)),1)/T;
    p_fly2(i,2) = p_fly1(bat_pairs(i,1))*p_fly1(bat_pairs(i,2));
end

%=== Calculate the difference in the number of flights across the session
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.2 0.9]);
tiledlayout(1,1,'TileSpacing','tight');
cum_flights = [zeros(1,n_tags); cumsum(max(diff(bflying),0))];
fligh_diff = zeros(T,n_pairs);
nexttile;
for i = 1:n_pairs
    fligh_diff(:,i) = cum_flights(:,bat_pairs(i,1))-cum_flights(:,bat_pairs(i,2));
    plot(t,fligh_diff(:,i),'LineWidth',2);    hold on;
    text(t(T),fligh_diff(T,i),[bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
end
rectangle('Position',[0 -10 t(T) 20],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
rectangle('Position',[0 -20 t(T) 40],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
hold off;
xlabel('Time (s)'); ylabel('Difference in the number of flights');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Analysis of the accelerometer Signal

figure('units','normalized','outerposition',[0.5 0 0.1 1]);
tiledlayout(n_tags,1,'TileSpacing','none');
for i=1:n_tags
    %=== Calculate FFT and PSD of the accelerometer signal
    extended_flight = movmax(bflying(:,i),[1*Fs 1*Fs]);
    trace = a_abs(find(extended_flight),i);
    n = 2^nextpow2(numel(trace));
    Y = fft(trace,n);
    P = abs(Y/n).^2;
    f = Fs*(0:n/2)/n;
    PSD = P(1:n/2+1);
    sm_PSD = smoothdata(PSD,'movmedian',n/Fs);
    %=== Find the peak corresponding to wingbeat
    [~,lo_f] = min(abs(f-6));
    [~,hi_f] = min(abs(f-10));
    [~,loc] = findpeaks(log(sm_PSD(lo_f:hi_f)),n/Fs,'MinPeakProminence',0.1);
    wb_F(i) = 6+max(loc);
    %=== Plot the spectrum and wingbeat harmonics
    nexttile;
    plot(f,PSD,'.k','LineWidth',2);     xlim([1 50]);    set(gca, 'YScale', 'log');
    hold on;
    plot(f,sm_PSD,'LineWidth',4,'Color',bat_clr(i,:));
    for j = 1:5,plot(wb_F(i)*[j j],ylim,'--k');end
    hold off;   xlabel('Hz');   yticks([]);
end
sgtitle('Accelerometer Spectrum');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Assign features associated with each bat and each pair

for i=1:n_tags   
    %=== Number of flights,p-value and percentage of clustered points, flight probability
    bat(i).f_num = f_num(i);
    bat(i).p_val_clus = p_val_clus(i,:);
    bat(i).CH_idx = CH_idx(i,:);
    bat(i).perc_clus = perc_clus(i,:);
    bat(i).p_fly1 = p_fly1(:,i);
    bat(i).wb_F = wb_F(:,i);
    
    %=== Exploration ratio see Obenhaus et al., "Functional network topography of the medial entorhinal cortex." PNAS 2022)
    occupancy = histcounts2(r(~bflying(:,i)& b_corr(:,i),1,i),r(~bflying(:,i)& b_corr(:,i),2,i),cell2mat(edges_d(1,1)),cell2mat(edges_d(1,2)));
    perimeter = [occupancy(1,:),occupancy(end,:),occupancy(:,1)',occupancy(:,end)'];
    bat(i).exp_ratio_all = nnz(occupancy>Fs*5)/numel(occupancy);
    bat(i).exp_ratio_per = nnz(perimeter>Fs*5)/numel(perimeter);
    bat(i).spatial_map = occupancy./T;

    %=== Fraction of time spent close to the feeder
    for j = 1:size(r_fd,1),condition(:,j) = vecnorm(r(:,:,i)-r_fd(j,:),2,2)<d_th;end
    bat(i).feedfraction = nnz(any(condition,2))/T;
    
    %=== Fraction of feeder,social,both and other flights
    bat(i).f_flights = (nnz(FLIGHTS.class(FLIGHTS.id==i)== 'f')+nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b'))/f_num(i);    % fraction of feeder flights, including social + feeder
    bat(i).s_flights = (nnz(FLIGHTS.class(FLIGHTS.id==i)== 's')+nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b'))/f_num(i);    % fraction of social flights, including social + feeder    
    bat(i).b_flights = nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b')/f_num(i);    
    bat(i).o_flights = nnz(FLIGHTS.class(FLIGHTS.id==i)== 'o')/f_num(i);
end

for i=1:n_pairs
    %=== Proximity/Landing index and p value
    pair(i).PI = PI(1,i);
    pair(i).PI_pval = PI_p(i);
    pair(i).CI = CI(1,i);
    pair(i).CI_pval = CI_p(i);
    pair(i).LI = LI(1,i,:);
    pair(i).LI_pval = LI_p(i,:);
    
    %=== Probability of simultaneous flight (observed vs calculated)
    pair(i).p_fly2 = p_fly2(i,:);
    
    %=== Maximal and median difference in the number of flights
    pair(i).max_flight_diff = max(abs(fligh_diff(:,i)));
    pair(i).med_flight_diff = median(fligh_diff(:,i));
    pair(i).sig_flight_diff = nnz(fligh_diff(:,i)>0)/T;
end

%% Save data

if options.save_data
    save([analysis_directory,'/Analyzed_Behavior_', batdate, '.mat'],...
            'bat','centroid','centroid_all','d_th','f_bd','f_clus','f_nn','f_smp','FLIGHTS','Group_name',...
            'n_ds','n_min','n_rep','options','pair','perc_time','pref_location','pref_time_loc','r_fd','samplesIn','T');           
end

%% Plot additional figures

%=== Show Density Histograms heat-map
figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
for i=1:n_tags
    sbpt = subplot(1,n_tags,i);
    hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
    xlabel('x');
    xlim(r_lim(1,:)); ylim(r_lim(2,:));   title(bat_nms(i,:));  view(90,90);  colormap(sbpt,custom_map(:,:,i)); % Change color scheme
    axis square;
end
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Show flights and positions during stationary epochs (cut the first 10s and last 3s)
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.2 0.35 0.3]);
tiledlayout(2,n_tags,'TileSpacing','tight');
for i = 1:n_tags
    nexttile(i);
    plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  view(90,90);
    nexttile(n_tags+i);
    scatter3(r_qt(10*Fs:end-3*Fs,1,i),r_qt(10*Fs:end-3*Fs,2,i),r_qt(10*Fs:end-3*Fs,3,i),5,'filled','MarkerFaceColor', bat_clr(i,:));
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
end
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Show all the clusters
figure;     
str = string(1:size(centroid_all,1));
textscatter3(centroid_all(:,1),centroid_all(:,2),centroid_all(:,3),str);    hold on;
scatter3(r_fd(:,1),r_fd(:,2),r_fd(:,3),50,'filled','MarkerFaceColor', 'r'); hold off;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Plot all centroids, with dimensions proportional to occupancy
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.2 0.3]);
for i = 1:n_tags
    scatter3(centroid{i,1}(:,1),centroid{i,1}(:,2),centroid{i,1}(:,3),samplesIn{i,1}./T*200,'filled');
    hold on;
end
axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   grid on;    hold off;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Plot Probabilities of flying together
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.2 0.33]);
scatter(p_fly2(:,2),p_fly2(:,1),'k','filled');
xlabel('P coupled fly (independent)');  ylabel('P coupled fly (real)');
text(p_fly2(:,2),p_fly2(:,1),bat_pair_nms,'VerticalAlignment','bottom','HorizontalAlignment','right')
axis square;    xlim(ylim);     hold on;    h = refline(1,0);   hold off;   h.Color = 'k';  h.LineStyle = '--';
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIG: Proximity index histograms
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,1,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(PI(:,i),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(PI(:))]);
    y1=get(gca,'ylim');  hold on; plot([PI(1,i) PI(1,i)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(PI_p(i),3)]);
end
xlabel('Proximity Index');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIG: Coupling index histograms
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,1,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(CI(:,i),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(CI(:))]);
    y1=get(gca,'ylim');  hold on; plot([CI(1,i) CI(1,i)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(CI_p(i),3)]);
end
xlabel('Coupling Index');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIG: Landing index histograms
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,2,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(LI(:,i,1),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(LI(:))]);
    y1=get(gca,'ylim');  hold on; plot([LI(1,i,1) LI(1,i,1)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(LI_p(i,1),3)]);
    ax(i+1) = nexttile;     histogram(LI(:,i,2),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(LI(:))]);
    y1=get(gca,'ylim');  hold on; plot([LI(1,i,2) LI(1,i,2)],y1); hold off;
    title([bat_nms(bat_pairs(i,2),:) '-' bat_nms(bat_pairs(i,1),:) ': ' num2str(LI_p(i,2),3)]);
end
xlabel('Proximity Index');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Network graph with PI
A = zeros(n_tags);
for i = 1:n_pairs
    if PI_p(i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif PI_p(i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif PI_p(i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;    
    end
end
G = graph(A,cellstr(bat_nms),'upper');
figure();
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); title('PI Network'); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Network graph with CI
A = zeros(n_tags);
for i = 1:n_pairs
    if CI_p(i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif CI_p(i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif CI_p(i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;    
    end
end
G = graph(A,cellstr(bat_nms),'upper');
figure();
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); title('CI Network'); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Network graph with LI
A = zeros(n_tags);
for i = 1:n_pairs
    if LI_p(i,1)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif LI_p(i,1)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif LI_p(i,1)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                  A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;    
    end
    if LI_p(i,2)<1e-3,     A(bat_pairs(i,2),bat_pairs(i,1)) = 0.7;
    elseif LI_p(i,2)<1e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.5;
    elseif LI_p(i,2)<5e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.3;
    else,                  A(bat_pairs(i,2),bat_pairs(i,1)) = 0.01;    
    end
end
G = digraph(A,cellstr(bat_nms));
figure();
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); title('LI Network'); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'ArrowSize',10,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Coordinated flight
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.2 0.3]);
for i = 1:n_tags
    plot(t(f_smp{i,1})/60,1:f_num(i),'Color', bat_clr(i,:),'LineWidth',3); hold on;
end
hold off;   xlabel('Time into session (min)');  ylabel('Flight number');    axis square;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

