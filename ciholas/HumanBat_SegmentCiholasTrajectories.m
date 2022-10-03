function [t,r_] = HumanBat_SegmentCiholasTrajectories(ciholas_data,ciholas2cortex,cortex_data,exp_data_path,batdate)

%% Function for extracting ciholas data and flights, based on A.F. script for the pre-processing of collective behavior (March 2022)

disp('...Processing Ciholas Data...');

%=== Parameters and metadata
n_tags = 1;
Fs = cortex_data.VideoFrameRate;                                                                                                      %Sampling frequency (Hz) for common time
bat_clr = lines(n_tags);                                                                                        %Bat Colors
v_th = 0.5;                                                                                                     %Velocity threshold (m/s) for flight segmentation

%=== Custom graded colormap(level,RGB,bat)
for i = 1:n_tags; for j = 1:3; custom_map(:,j,i) = linspace(1,bat_clr(i,j))'; end; end

%% Remove duplicate samples and shift Ciholas Time by 1.05s (t = 0 corresponds to the first '3s' Master-9 TTL)

[~,ia,~] = unique(ciholas_data.tag_data{7}(:,8),'stable');         ciholas_data.tag_data{7} = ciholas_data.tag_data{7}(ia,:);              ciholas_data.tag_data{7}(:,8) = ciholas_data.tag_data{7}(:,8)+1.05;
[~,ia,~] = unique(ciholas_data.tag_data_filt{7}(:,8),'stable');    ciholas_data.tag_data_filt{7} = ciholas_data.tag_data_filt{7}(ia,:);    ciholas_data.tag_data_filt{7}(:,8) = ciholas_data.tag_data_filt{7}(:,8)+1.05;
[~,ia,~] = unique(ciholas_data.tag_ac_data{7}(:,8),'stable');      ciholas_data.tag_ac_data{7} = ciholas_data.tag_ac_data{7}(ia,:);        ciholas_data.tag_ac_data{7}(:,8) = ciholas_data.tag_ac_data{7}(:,8)+1.05;

%% Calculate evenly sampled kinematic variables (r,v,a) and wing-beats epochs

% Find length of cortex TTl vector
cortex_ttls = cortex_data.AnalogSignals(:,2);
[cortex_ttl_val,cortex_ttl_idx] = findpeaks(cortex_ttls,'MinPeakHeight',4.5,'MinPeakDistance',300);
first_cortex_ttl = cortex_ttl_idx(1);       last_cortex_ttl = cortex_ttl_idx(end);
length_of_vector_to_cover_all_cortex_ttls = length(cortex_ttls(first_cortex_ttl:last_cortex_ttl))/Fs;

% t - a vector the length of the number of samples covered by the cortex TTLs at 120Hz
t = [0:1/Fs:length_of_vector_to_cover_all_cortex_ttls]';       %Evenly sampled time vector from first to last TTL.
T = length(t);                                                     %Number of time samples
r = zeros(T,3,n_tags);                                             %3D position vector (sample,dim,id)

%=== Interpolate position at evenly spaced time points sampled to 120Hz
r =  interp1(ciholas_data.tag_data_filt{7}(:,8), ciholas_data.tag_data_filt{7}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!

%=== Interpolate acceleration at evenly spaced time points
a = zeros(T,3,n_tags);
a(:,:) =  interp1(ciholas_data.tag_ac_data{7}(:,8), ciholas_data.tag_ac_data{7}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
a_abs = squeeze(vecnorm(a,2,2));   %Modulus
a_flt = bandpass(a_abs,[7 9],Fs);  %Filtered at the wing-beat frequency

%=== Transform ciholas position and acceleration to cortex coordinates
r_ = ciholas2cortex.*r;             
a_ = ciholas2cortex.*a;

%=== Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v,2,2));  angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%=== Detect flight epochs based on wing-beat signal
figure(); hold on;
[up,lo] = envelope(a_flt+normrnd(0,1e-3,length(a_flt),1),Fs/10,'peak');   %Envelope of the acceleration signal (noise addedd to avoid problems with splining)
env = normalize(up - lo,'range');                                           %Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                                       %Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,1.67*Fs)>Fs/6; 
%figure(); hold on; %Heuristic criterion for flight detection
     area(t,wBeats*3,'FaceAlpha',0.3,'LineStyle','none');  hold on;
     plot(t,normalize(v_abs,'range'));
     plot(t,r(:,1),t,r(:,2));
     plot(t,normalize(a_flt,'range',[-1 1]));
     plot(t,normalize(movsum(env>env_th,1.67*Fs),'range'));
     plot(t,up);
     plot(t,lo);
     hold off;

% % Uncomment to control interpolation quality
% ax(1) = subplot(311);  plot(ciholas_data.tag_data{1,4}(:,8),ciholas_data.tag_data{1,4}(:,3)','.',t,r(:,1,4),'-');
% ax(2) = subplot(312);  plot(ciholas_data.tag_data{1,4}(:,8),ciholas_data.tag_data{1,4}(:,4)','.',t,r(:,2,4),'-');
% ax(3) = subplot(313);  plot(ciholas_data.tag_data{1,4}(:,8),ciholas_data.tag_data{1,4}(:,5)','.',t,r(:,3,4),'-');
% linkaxes(ax,'x');
% ax(1) = subplot(311);  plot(ciholas_data.tag_ac_data{1,1}(:,8),ciholas_data.tag_ac_data{1,1}(:,3)','.',t,a(:,1,1),'-');
% ax(2) = subplot(312);  plot(ciholas_data.tag_ac_data{1,1}(:,8),ciholas_data.tag_ac_data{1,1}(:,4)','.',t,a(:,2,1),'-');
% ax(3) = subplot(313);  plot(ciholas_data.tag_ac_data{1,1}(:,8),ciholas_data.tag_ac_data{1,1}(:,5)','.',t,a(:,3,1),'-');
% linkaxes(ax,'x');

%% Correct position by median filtering when the bat is not flapping its wings

    %=== Find stationary epochs
    stat_periods = repmat(~wBeats,1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_stat = r;    r_stat(~stat_periods) = nan;
    
    %=== Filter position during stationary epochs
    r_stat =  smoothdata(r_stat(:,:),1,'movmedian',Fs*6,'omitnan');
    r_stat(~stat_periods) = nan;
    
    % % Uncomment to control filtering quality
%      ax(1) = subplot(311);  plot(tag_data{7}(:,8),tag_data{7}(:,3)',t,r_stat(:,1,7));
%      ax(2) = subplot(312);  plot(tag_data{7}(:,8),tag_data{7}(:,4)',t,r_stat(:,2,7));
%      ax(3) = subplot(313);  plot(tag_data{7}(:,8),tag_data{7}(:,5)',t,r_stat(:,3,7));
%      linkaxes(ax,'x');
    
    %=== Substitute median filtered data when bat is not flapping its wings
    r_stat = fillmissing(r_stat,'constant',0,1);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);     %DO NOT USE lowess!!
    
%      ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)',t,r_corr(:,1,4));
%      ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)',t,r_corr(:,2,4));
%      ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)',t,r_corr(:,3,4));
%      linkaxes(ax,'x');
    
% else
%     %Aggressive median filtering on original data
%     tag_data_stat = tag_data;
%     r_stat = zeros(T,3,n_tags);
%     for i = 1:n_tags
%         tag_data_stat{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'movmedian',Fs*3);
%         r_stat(:,:,i) =  csaps(tag_data_stat{1,i}(:,8), tag_data_stat{1,i}(:,[3:5])', 1, t)';
%     end
%     %Calculate xy-velocity and smooth
%     v_filt = smoothdata(squeeze( vecnorm( v(:,1:2,:),2,2)), 1, 'movmedian', Fs*3);
%     %Substitute median filtered data when velocity is less than threshold
%     stat_periods = repmat((v_filt < v_th),1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
%     r_corr = r_stat.*stat_periods+r.*(~stat_periods);
%     r_corr = smoothdata(r_corr,1,'loess',Fs*1);
%     r_corr_1 = smoothdata(r_corr,1,'loess',Fs*1);
% end

%% Correct position, recalculate velocity and heading angle

% if options.use_r_corr
%     r_old = r;  v_old = v_abs;  %Store old versions, in case you need them
%     r = r_corr; 
%     v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v(:,1:3,:),2,2)); angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));
% end

%% Flight segmentation

%=== Detect flight starts and stops by using risetime and falltime
%!!! Be careful with the levels of risetime and falltime, these
%    are fundamental in order to exclude stationary tails in the flight
bflying = zeros(T,n_tags); f_num = zeros(n_tags,1);   f_smp = cell(n_tags,2);
[bflying,f_num,f_smp{1},f_smp{2}] = FlightSegm_AF_v0(v_abs,v_th,Fs);

%Define staionary periods when the bat is not flying
stat_periods = repmat(~bflying,1,1,3);   stat_periods = permute(stat_periods,[1 3 2]);
r_qt = r;   r_qt(~stat_periods)=nan;     angle(squeeze(stat_periods(:,1,:))) = nan;

%% Make a few figures

    
%=== FIGURE: Raw Position VS Corrected Position
for j = 1:3
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(ciholas_data.tag_data{7}(:,8), ciholas_data.tag_data{7}(:,2+j),t,r_corr(:,j));
    sgtitle(['D-' num2str(j)]);   legend('raw','corrected');    ylabel('m');     xlabel('Time(s)'); xlim([0,t(T)]);
end
 
%=== FIGURE: Raw Velocity VS Corrected Velocity
% figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% plot(t,v_old,'.');    hold on;     sgtitle('v');
% plot(t,v_abs,'.');    hold off;    legend('raw','corrected');      ylabel('m/s');   xlabel('Time(s)'); xlim([0,t(T)]);  

%=== FIGURE: Scatter plot all bats
figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
subplot(131);  plot3(r(:,1),r(:,2),r(:,3),'Color', bat_clr);  title('3D view');                 hold on;  axis equal;
subplot(132);  plot3(r(:,1),r(:,2),r(:,3),'Color', bat_clr);   title('Top view');   view(0,90);  hold on;  axis equal;
subplot(133);  plot3(r(:,1),r(:,2),r(:,3),'Color', bat_clr);   title('Door view');  view(-30,0); hold on;  axis equal;
hold off;

%=== FIGURE: Density Histograms heat-map
figure(); 
hist3(r(:,1:2),'CdataMode','auto','edgecolor','none','FaceColor','interp');
xlabel('x');
view(90,90);  

%=== Polar plots of heading direction 


%=== FIGURE: Velocity and flight segmentation
figure();       set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(n_tags,1,'TileSpacing','tight');
area(t,bflying*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
area(t,wBeats*-1,'FaceAlpha',0.3,'LineStyle','none');  refline(0,-0.3);
plot(t,v_abs,'.','Color', bat_clr);     plot(t,r(:,1),'k--');  ylabel('Velocity (m/s)');     hold off;
legend('Fly','Wing-B','Vel','x(m)');
xlabel('Time(s)');

%=== Pass back a struct that contains the 3D positions in Cortex 120 in
%Cortex coordinates of the flights from ciholas
% Save data
save([strcat(exp_data_path,'ciholas/'),'/Extracted_Behavior_', num2str(batdate), '.mat'],...
    'a','a_abs','a_flt','angle','bat_clr',...
    'batdate','bflying','env',...
    'f_num','f_smp','Fs','n_tags',...
    'r','r_qt','stat_periods',...
    't','T','v','v_abs','v_th','wBeats');


%% Save figures and data
% if options.save_data
%     figHandles = findall(0,'Type','figure');
%     for i = 1:numel(figHandles)
%         saveas(figHandles(i),[analysis_directory, '/', batdate '_figure' num2str(numel(figHandles)+1-i) '.png']);
%     end
%     close all;
%     save([analysis_directory,'/Extracted_Behavior_', batdate, '.mat'],...
%         'a','a_abs','a_flt','angle','bat_clr','bat_nms','bat_pair_nms','bat_pairs',...
%         'batdate','bflying','CDPmtdata','edges_d','env','extracted_CDPfile',...
%         'f_num','f_smp','Fs','Group_name','n_tags','options',...
%         'r','r_lim','r_old','r_qt','stat_periods',...
%         't','T','v','v_abs','v_th','wBeats');
% end
% 
% %% Save movie at 10Hz with bat trajectories
% if options.savemovie
%     saveCollMovie_AF_v0(r,Fs,bat_nms,r_lim,analysis_directory,batdate);
% end
end

