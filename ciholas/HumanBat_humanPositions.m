function [posK,posM,t,r_,M_flights,K_flights] = HumanBat_humanPositions(ciholas_data,ciholas2cortex,cortex_data,exp_data_path,batdate)

%HumanBat_humanPositions Summary of this function goes here
%   Parameters
%   ----------
%   extracted_ciholas_mat
%       Path to directory that contains ciholas/ processed data
%
%
%   Outputs
%   -------
%   walkpaths 
%       struct of human positions

% kq_rh = ciholas_data.tag_data_filt{1}(:,3:5);
% ms_rh = ciholas_data.tag_data_filt{2}(:,3:5);
% kq_rp = ciholas_data.tag_data_filt{3}(:,3:5);
% ms_rp = ciholas_data.tag_data_filt{4}(:,3:5);
% ms_lh = ciholas_data.tag_data_filt{5}(:,3:5);
% kq_lh = ciholas_data.tag_data_filt{6}(:,3:5);
% blondie = ciholas_data.tag_data_filt{7}(:,3:5);
% seconds = ciholas_data.tag_data_filt{7}(:,8);

% Ground truth of tripods
tripods = [-2.4,0.56,1.1;
           -1.6,2.07,1.1;
           2.07,1.6,0.9;
           2.26,0.1,0.9;
           2.1,-1.6,1;
           0.36,-2.2,0.9];
ntags = 6;
Fs = cortex_data.VideoFrameRate;                                                                                                      %Sampling frequency (Hz) for common time
bat_clr = lines(ntags);                                                                                        %Bat Colors
v_th = 0.5; 

%% Resample the ciholas data
% Remove duplicate samples and shift Ciholas Time by 1.05s (t = 0 corresponds to the first '3s' Master-9 TTL)
for i=1:ntags
    [~,ia,~] = unique(ciholas_data.tag_data{1,i}(:,8),'stable');         ciholas_data.tag_data{1,i} = ciholas_data.tag_data{1,i}(ia,:);              ciholas_data.tag_data{1,i}(:,8) = ciholas_data.tag_data{1,i}(:,8)+1.05;
    [~,ia,~] = unique(ciholas_data.tag_data_filt{1,i}(:,8),'stable');    ciholas_data.tag_data_filt{1,i} = ciholas_data.tag_data_filt{1,i}(ia,:);    ciholas_data.tag_data_filt{1,i}(:,8) = ciholas_data.tag_data_filt{1,i}(:,8)+1.05;
    [~,ia,~] = unique(ciholas_data.tag_ac_data{1,i}(:,8),'stable');      ciholas_data.tag_ac_data{1,i} = ciholas_data.tag_ac_data{1,i}(ia,:);        ciholas_data.tag_ac_data{1,i}(:,8) = ciholas_data.tag_ac_data{1,i}(:,8)+1.05;
end
%% Calculate evenly sampled kinematic variables (r,v,a) and wing-beats epochs

% Find length of cortex TTl vector
cortex_ttls = cortex_data.AnalogSignals(:,2);
[cortex_ttl_val,cortex_ttl_idx] = findpeaks(cortex_ttls,'MinPeakHeight',4.5,'MinPeakDistance',300);
first_cortex_ttl = cortex_ttl_idx(1);       last_cortex_ttl = cortex_ttl_idx(end);
length_of_vector_to_cover_all_cortex_ttls = length(cortex_ttls(first_cortex_ttl:last_cortex_ttl))/Fs;

% t - a vector the length of the number of samples covered by the cortex TTLs at 120Hz
t = [0:1/Fs:length_of_vector_to_cover_all_cortex_ttls]';       %Evenly sampled time vector from first to last TTL. Seconds 
T = length(t);                                                     %Number of time samples
r = zeros(T,3,ntags);                                             %3D position vector (sample,dim,id)

%=== Interpolate position at evenly spaced time points (already resampled)
for i = 1:ntags
    r(:,:,i) =  interp1(ciholas_data.tag_data_filt{1,i}(:,8), ciholas_data.tag_data_filt{1,i}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
end

%=== Interpolate acceleration at evenly spaced time points
a = zeros(T,3,ntags);
for i = 1:ntags
    a(:,:,i) =  interp1(ciholas_data.tag_ac_data{1,i}(:,8), ciholas_data.tag_ac_data{1,i}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
end
a_abs = squeeze(vecnorm(a,2,2));    %Modulus
a_flt = bandpass(a_abs,[7 9],Fs);  %Filtered at the wing-beat frequency

%=== Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,ntags),v);   v_abs = squeeze(vecnorm(v,2,2));  angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%=== Transform ciholas position and acceleration to cortex coordinates
r_ = ciholas2cortex.*r;             
a_ = ciholas2cortex.*a;

%=== Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,ntags),v);   v_abs = squeeze(vecnorm(v,2,2));  angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%% Perform FFT and Bessel filter on the human data to find positions

% 1. Take fft of small segment of raw data to find frequencies, build
% filter that filters out BELOW those frequencies (If not in cycles, multiply by 2pi/sample rate)
figure(); plot(r_(:,1,2));
x_seg = r_(6000:400000,1,2);      y_seg = r_(6000:400000,2,2);       z_seg = r_(6000:400000,3,2);  
%x_seg = r_(74000:99000,1,2);ciholas_data.tag_data_filt{2}(74000:99000,3);      y_seg = ciholas_data.tag_data_filt{2}(74000:99000,4);      z_seg = ciholas_data.tag_data_filt{2}(74000:99000,5);

% 2a. Normalize raw data
x_seg_z = normalize(x_seg);     y_seg_z = normalize(y_seg);     z_seg_z = normalize(z_seg); 

% 2b. Take FFT
xdft_x = fft(x_seg_z);    xdft_x = xdft_x(1:length(x_seg_z)/2+1);     psdx_x = (1/(100*length(x_seg_z)))*abs(xdft_x).^2;      psdx_x(2:end-1) = 2*psdx_x(2:end-1);     xdft_y = fft(y_seg_z);    xdft_y = xdft_y(1:length(y_seg_z)/2+1);     psdx_y = (1/(100*length(y_seg_z)))*abs(xdft_y).^2;      psdx_y(2:end-1) = 2*psdx_y(2:end-1);   xdft_z = fft(z_seg_z);    xdft_z = xdft_z(1:length(z_seg_z)/2+1);     psdx_z = (1/(100*length(z_seg_z)))*abs(xdft_z).^2;      psdx_z(2:end-1) = 2*psdx_z(2:end-1);
freq_x = 0:100/length(x_seg_z):100/2;    freq_y = 0:100/length(y_seg_z):100/2;   freq_z = 0:100/length(z_seg_z):100/2;

% 3. Besel (increase order if artifacts) 
cf = 1.3;    [z,p,k] = besself(4,cf,'low');      [zd,pd,kd] = bilinear(z,p,k,100);     [sos,g] = zp2sos(zd,pd,kd);

% 4. Filt Filt (forces linear phase)
x_seg_filtfilt = filtfilt(sos,g,x_seg_z);   y_seg_filtfilt = filtfilt(sos,g,y_seg_z);   z_seg_filtfilt = filtfilt(sos,g,z_seg_z);

% Plot the raw, fft, and filtered segment of data
figure(); hold on; title("Small segment of jitter when human #1 is still"); 
subplot(3,1,1); plot(x_seg_z);
subplot(3,1,2); plot(freq_x,10*log10(psdx_x)); %plot(f,P1); 
subplot(3,1,3); plot(x_seg_filtfilt);

% === Apply this to the whole data
htag = [1,2,3,4,5,6];
for i=1:length(htag)
    D_raw_x{i} = r_(:,1,htag(i)); D_raw_y{i} = r_(:,2,htag(i));  D_raw_z{i} = r_(:,3,htag(i)); 
    %D_raw_x{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,3); D_raw_y{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,4); D_raw_z{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,5);

    % normalize 
    D_z_x{i} = normalize(D_raw_x{i}); D_z_y{i} = normalize(D_raw_y{i}); D_z_z{i} = normalize(D_raw_z{i});
    
    % fft
    xdft_x = fft(D_z_x{i});    xdft_x = xdft_x(1:length(D_z_x{i})/2+1);     psdx_x = (1/(100*length(D_z_x{i})))*abs(xdft_x).^2;      psdx_x(2:end-1) = 2*psdx_x(2:end-1);   xdft_y = fft(D_z_y{i});    xdft_y = xdft_y(1:length(D_z_y{i})/2+1);     psdx_y = (1/(100*length(D_z_y{i})))*abs(xdft_y).^2;      psdx_y(2:end-1) = 2*psdx_y(2:end-1);     xdft_z = fft(D_z_z{i});    xdft_z = xdft_z(1:length(D_z_z{i})/2+1);     psdx_z = (1/(100*length(D_z_z{i})))*abs(xdft_z).^2;      psdx_z(2:end-1) = 2*psdx_z(2:end-1);
    freq_x = 0:100/length(D_z_x{i}):100/2;   freq_y = 0:100/length(D_z_y{i}):100/2;   freq_z = 0:100/length(D_z_z{i}):100/2;
    
    % filter with bessel
    D_z_x_filtfilt{i} = filtfilt(sos,g,D_z_x{i});    D_z_y_filtfilt{i} = filtfilt(sos,g,D_z_y{i});    D_z_z_filtfilt{i} = filtfilt(sos,g,D_z_z{i}); 

    % filter with movmed
    D_z_x_filtfilt_movmed{i} = movmedian(D_z_x_filtfilt{i},[4000 4000]);     D_z_y_filtfilt_movmed{i} = movmedian(D_z_y_filtfilt{i},[4000 4000]);     D_z_z_filtfilt_movmed{i} = movmedian(D_z_z_filtfilt{i},[4000 4000]); 
    
    % plot
    figure(); hold on; title("Filtered vs Raw Data Ciholas Human #1"); 
    subplot(4,1,1); plot(D_raw_x{i});
    subplot(4,1,2); plot(freq_x,10*log10(psdx_x)); %plot(f,P1); 
    subplot(4,1,3); plot(D_z_x_filtfilt{i});
    subplot(4,1,4); hold on; plot(D_z_x{i}); plot(D_z_x_filtfilt{i}); hold off;
    hold off;
end

% Try applying a movmed filter too
for i=1:length(htag)
    figure(); hold on;
    %plot(D_raw{i}); plot(D_z_filtfilt{i}); plot(D_z_filtfilt_movmed{i}); hold off;
    plot(diff(D_z_x_filtfilt_movmed{i})*1300);
    plot(D_z_x_filtfilt_movmed{i});
end

% === Assign to posK and posM
posK = [D_raw_x{3},D_raw_y{3},D_raw_z{3}];
posM = [D_raw_x{4},D_raw_y{4},D_raw_z{4}];

% === Extract static positions
% Split ciholas acceleration data to find human positions during the first
% 60 minutes

% human_locs contains the aligned timestamps at 120Hz of the
% beginnings/endings of the human flights 
human_locs = {}; D_z_x_filtfilt_movmed_diff = {}; D_z_y_filtfilt_movmed_diff = {}; 
figure(); hold on; 
for i=1:length(htag)
    clear pk_vect yy_pks yy_locs yy
    % take derivative
    D_z_x_filtfilt_movmed_diff{i} = diff(D_z_x_filtfilt_movmed{i})*1000;
    D_z_y_filtfilt_movmed_diff{i} = diff(D_z_y_filtfilt_movmed{i})*1000;
    [yy_pks,yy_locs] = findpeaks(abs(D_z_y_filtfilt_movmed_diff{i}),'MinPeakDistance',10000,'MinPeakProminence',1);
    pk_vect = NaN(1,length(abs(D_raw_x{i})));
    pk_vect(yy_locs) = 2500;
    subplot(2,3,i); hold on; title(strcat("Human #",num2str(i)," ","switches.")); plot(D_z_x_filtfilt_movmed{i}); plot(D_raw_x{i}); stem(pk_vect); hold off;
    human_locs{i} = yy_locs;
    linkaxes;
end

% % === Extract human flights 

% Find the envalope of the y trace
[y_up,y_lo] = envelope(abs(D_z_y_filtfilt_movmed_diff{i}),2000,'peak');

[yy_e_pks,yy_e_locs] = findpeaks(y_up,'MinPeakProminence',1.4);
pk_vect_e = NaN(1,length(y_up));
pk_vect_e(yy_e_locs) = 2500;
figure(); hold on; plot(D_z_y_filtfilt_movmed{i}); plot(D_raw_y{i}); stem(pk_vect_e); hold off;
figure(); hold on; plot(abs(D_z_y_filtfilt_movmed{i})); plot(y_up);

%for i=1:length(yy_e_locs)
%    zero_range_raw = D_raw_y{i}(yy_e_locs(i)-1000:yy_e_locs(i)+2000);
%    zero_range_env = D_z_y_filtfilt_movmed{i}(yy_e_locs(i)-1000:yy_e_locs(i)+2000);
%    figure(); hold on; plot(D_z_y_filtfilt_movmed{i}); scatter(yy_e_locs(i)-1000:yy_e_locs(i)+2000,zero_range_env);
%end

%% MADELEINE: For each peak found, plot the "flight" until the next peak

% madeleine_tag = 2;
% extra = 1;
% for i=1:length(yy_locs)-1
%     figure(); 
%     subplot(2,1,1); hold on;
%     plot3(D_raw_x{madeleine_tag},D_raw_y{madeleine_tag},D_raw_z{madeleine_tag},'Color',[0.8 0.8 0.8]);
%     plot3(D_raw_x{madeleine_tag}(yy_locs(i):yy_locs(i+extra)),D_raw_y{madeleine_tag}(yy_locs(i):yy_locs(i+extra)),D_raw_z{madeleine_tag}(yy_locs(i):yy_locs(i+extra)));
%     xlim([-2.6*1e3 2.6*1e3]); ylim([-2.6*1e3 2.6*1e3]);
%     scatter3(D_raw_x{madeleine_tag}(yy_locs(i)),D_raw_y{madeleine_tag}(yy_locs(i)),D_raw_z{madeleine_tag}(yy_locs(i)),200,'b*');
%     scatter3(D_raw_x{madeleine_tag}(yy_locs(i+extra)),D_raw_y{madeleine_tag}(yy_locs(i+extra)),D_raw_z{madeleine_tag}(yy_locs(i+extra)),200,'b*'); hold off;
%     subplot(2,1,2);
%     plot3(diff())
% end

% === For each human, for each position, plot the ciholas tags from that
% position to make sure it's all around the same area
% m = 1*Fs*60;
% for i=1:2
%     for j=1:length(human_locs{i})
%         % Get the median position during that switch
%         human_med_positions(j,:,i) = [median(ciholas_data.tag_data_filt{htag(i)}(human_locs{i}(j):human_locs{i}(j)+m,3)),...
%                                     median(ciholas_data.tag_data_filt{htag(i)}(human_locs{i}(j):human_locs{i}(j)+m,4)),...
%                                     median(ciholas_data.tag_data_filt{htag(i)}(human_locs{i}(j):human_locs{i}(j)+m,5))];
%     end
%     figure(); hold on; axis([-2.6 2.6 -2.6 2.6 0 Inf]); plot3(ciholas_data.tag_data_filt{htag(i)}(human_locs{i}(j):human_locs{i}(j)+m,3),ciholas_data.tag_data_filt{htag(i)}(human_locs{i}(j):human_locs{i}(j)+m,4),ciholas_data.tag_data_filt{htag(i)}(human_locs{i}(j):human_locs{i}(j)+m,5)); hold off;
% end
% hold off;

% == For each human, plot the order of tripods visited
% human_tripod_vector = [];
% for i=1:2
%     for j=1:length(D_z_x_filtfilt_movmed{i})
%         if (tripods(1,1)+0.6 < D_z_x_filtfilt_movmed{i}(j) < tripods(1,1)-0.6 && 2==2)
%             human_tripod_vector(j) = 1;
%         elseif (tripods(2,1)+0.6 < D_z_x_filtfilt_movmed{i}(j) < tripods(2,1)-0.6 && 2==2)
%             human_tripod_vector(j) = 1;
%         end
%     end
% end

M_flights = [];
K_flights = [];
        

% figure(); hold on; title("Average Human positions (red Kevin blue Madeleine)"); scatter3(human_med_positions(:,1,1),human_med_positions(:,2,1),human_med_positions(:,3,1),'r');
% scatter3(human_med_positions(1:9,1,2),human_med_positions(1:9,2,2),human_med_positions(1:9,3,2),'b');