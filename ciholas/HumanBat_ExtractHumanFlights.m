function [Mstart,Mend,Kstart,Kend] = HumanBat_ExtractHumanFlights(human_r,human_t)

%=== Calculate velocity and 2D-direction of motion for Madeleine
human_M = squeeze(human_r(:,:,4)); human_K = squeeze(human_r(:,:,3));

v = diff(human_M,1,1)./diff(human_t); v=cat(1,zeros(1,3,1),v);   v_abs = squeeze(vecnorm(v,2,2));  angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%% Perform FFT and Bessel filter on the human data to find positions

% 1. Take fft of small segment of raw data to find frequencies, build
% filter that filters out BELOW those frequencies (If not in cycles, multiply by 2pi/sample rate)
x_seg = human_M(6000:400000,1);      y_seg = human_M(6000:400000,2);       z_seg = human_M(6000:400000,3);  

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

%% === Apply this to Madeleine
D_raw_x = human_M(:,1); D_raw_y = human_M(:,2);  D_raw_z = human_M(:,3); 
%D_raw_x{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,3); D_raw_y{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,4); D_raw_z{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,5);

% normalize 
D_z_x = normalize(D_raw_x); D_z_y = normalize(D_raw_y); D_z_z = normalize(D_raw_z);

% fft
xdft_x = fft(D_z_x);    xdft_x = xdft_x(1:length(D_z_x)/2+1);     psdx_x = (1/(100*length(D_z_x)))*abs(xdft_x).^2;      psdx_x(2:end-1) = 2*psdx_x(2:end-1);   xdft_y = fft(D_z_y);    xdft_y = xdft_y(1:length(D_z_y)/2+1);     psdx_y = (1/(100*length(D_z_y)))*abs(xdft_y).^2;      psdx_y(2:end-1) = 2*psdx_y(2:end-1);     xdft_z = fft(D_z_z);    xdft_z = xdft_z(1:length(D_z_z)/2+1);     psdx_z = (1/(100*length(D_z_z)))*abs(xdft_z).^2;      psdx_z(2:end-1) = 2*psdx_z(2:end-1);
freq_x = 0:100/length(D_z_x):100/2;   freq_y = 0:100/length(D_z_y):100/2;   freq_z = 0:100/length(D_z_z):100/2;

% filter with bessel
D_z_x_filtfilt = filtfilt(sos,g,D_z_x);    D_z_y_filtfilt = filtfilt(sos,g,D_z_y);    D_z_z_filtfilt = filtfilt(sos,g,D_z_z); 

% filter with movmed
D_z_x_filtfilt_movmed = movmedian(D_z_x_filtfilt,[4000 4000]);     D_z_y_filtfilt_movmed = movmedian(D_z_y_filtfilt,[4000 4000]);     D_z_z_filtfilt_movmed = movmedian(D_z_z_filtfilt,[4000 4000]); 

% plot
figure(); hold on; title("Filtered vs Raw Data Ciholas Human #1"); 
subplot(4,1,1); plot(D_raw_x);
subplot(4,1,2); plot(freq_x,10*log10(psdx_x)); %plot(f,P1); 
subplot(4,1,3); plot(D_z_x_filtfilt);
subplot(4,1,4); hold on; plot(D_z_x); plot(D_z_x_filtfilt); hold off;
hold off;

% Apply movmed filter too
figure(); hold on;
plot(diff(D_z_x_filtfilt_movmed)*1300);
plot(D_z_x_filtfilt_movmed);

% === Assign to posK and posM
posM = [D_raw_x,D_raw_y,D_raw_z];

% === Extract static positions
figure(); hold on; 
clear pk_vect yy_pks yy_locs yy
% take derivative
D_z_x_filtfilt_movmed_diff = diff(D_z_x_filtfilt_movmed)*1000;
D_z_y_filtfilt_movmed_diff = diff(D_z_y_filtfilt_movmed)*1000;
[yy_pks,yy_locs] = findpeaks(abs(D_z_y_filtfilt_movmed_diff),'MinPeakDistance',100,'MinPeakProminence',1);
pk_vect = NaN(1,length(abs(D_raw_x)));
pk_vect(yy_locs) = 2500;
title(strcat("Human M switches.")); plot(D_z_x_filtfilt_movmed); plot(D_raw_x); stem(pk_vect); hold off;
M_locs = yy_locs;
linkaxes;

%Start of stillness is pk_vect - 100 samples
Mp = find(~isnan(pk_vect));
for i=1:length(Mp)-1
    Mstart(i) = Mp(i)+600;
    Mend(i) = Mp(i+1)-1000;
end

for i=1:length(Mstart)-1
    tc=i;
    figure(); hold on; title(strcat("M ",num2str(tc))); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]); scatter3(human_r(Mstart(tc):Mend(tc),1,4),human_r(Mstart(tc):Mend(tc),2,4),human_r(Mstart(tc):Mend(tc),3,4));
end

% Clean manually ugh
end_clip = 0;
start_clip = 0;
tc=3
figure(); hold on; title(strcat("M ",num2str(tc))); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]); 
scatter3(human_r(Mstart(tc)+start_clip:Mend(tc)-end_clip,1,4),human_r(Mstart(tc)+start_clip:Mend(tc)-end_clip,2,4),human_r(Mstart(tc)+start_clip:Mend(tc)-end_clip,3,4));


%% === Apply this to Kevin
D_raw_x = human_K(:,1); D_raw_y = human_K(:,2);  D_raw_z = human_K(:,3); 
%D_raw_x{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,3); D_raw_y{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,4); D_raw_z{i} = ciholas_data.tag_data{htag(i)}(1:length(ciholas_data.tag_data{htag(i)})-30000,5);

% normalize 
D_z_x = normalize(D_raw_x); D_z_y = normalize(D_raw_y); D_z_z = normalize(D_raw_z);

% fft
xdft_x = fft(D_z_x);    xdft_x = xdft_x(1:length(D_z_x)/2+1);     psdx_x = (1/(100*length(D_z_x)))*abs(xdft_x).^2;      psdx_x(2:end-1) = 2*psdx_x(2:end-1);   xdft_y = fft(D_z_y);    xdft_y = xdft_y(1:length(D_z_y)/2+1);     psdx_y = (1/(100*length(D_z_y)))*abs(xdft_y).^2;      psdx_y(2:end-1) = 2*psdx_y(2:end-1);     xdft_z = fft(D_z_z);    xdft_z = xdft_z(1:length(D_z_z)/2+1);     psdx_z = (1/(100*length(D_z_z)))*abs(xdft_z).^2;      psdx_z(2:end-1) = 2*psdx_z(2:end-1);
freq_x = 0:100/length(D_z_x):100/2;   freq_y = 0:100/length(D_z_y):100/2;   freq_z = 0:100/length(D_z_z):100/2;

% filter with bessel
D_z_x_filtfilt = filtfilt(sos,g,D_z_x);    D_z_y_filtfilt = filtfilt(sos,g,D_z_y);    D_z_z_filtfilt = filtfilt(sos,g,D_z_z); 

% filter with movmed
D_z_x_filtfilt_movmed = movmedian(D_z_x_filtfilt,[4000 4000]);     D_z_y_filtfilt_movmed = movmedian(D_z_y_filtfilt,[4000 4000]);     D_z_z_filtfilt_movmed = movmedian(D_z_z_filtfilt,[4000 4000]); 

% plot
figure(); hold on; title("Filtered vs Raw Data Ciholas Human Kevin"); 
subplot(4,1,1); plot(D_raw_x);
subplot(4,1,2); plot(freq_x,10*log10(psdx_x)); %plot(f,P1); 
subplot(4,1,3); plot(D_z_x_filtfilt);
subplot(4,1,4); hold on; plot(D_z_x); plot(D_z_x_filtfilt); hold off;
hold off;

% Apply movmed filter too
figure(); hold on;
plot(diff(D_z_x_filtfilt_movmed)*1300);
plot(D_z_x_filtfilt_movmed);

% === Assign to posK and posM
posK = [D_raw_x,D_raw_y,D_raw_z];

% === Extract static positions

figure(); hold on; 
clear pk_vect yy_pks yy_locs yy
% take derivative
D_z_x_filtfilt_movmed_diff = diff(D_z_x_filtfilt_movmed)*1000;
D_z_y_filtfilt_movmed_diff = diff(D_z_y_filtfilt_movmed)*1000;
[yy_pks,yy_locs] = findpeaks(abs(D_z_y_filtfilt_movmed_diff),'MinPeakDistance',100,'MinPeakProminence',1);
pk_vect = NaN(1,length(abs(D_raw_x)));
pk_vect(yy_locs) = 2500;
title(strcat("Human K switches.")); plot(D_z_x_filtfilt_movmed); plot(D_raw_x); stem(pk_vect); hold off;
K_locs = yy_locs;
linkaxes;

%Start of stillness is pk_vect - 100 samples
Kp = find(~isnan(pk_vect));
for i=1:length(Kp)-1
    Kstart(i) = Kp(i)+600;
    Kend(i) = Kp(i+1)-1000;
end

% Display all individual perdios of stillness and manually list which ones
% to use
for i=1:length(Kstart)-1
    tc=i;
    figure(); hold on; title(strcat("K ",num2str(tc))); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]); scatter3(human_r(Kstart(tc):Kend(tc),1,3),human_r(Kstart(tc):Kend(tc),2,3),human_r(Kstart(tc):Kend(tc),3,3));
end

% Clean manually ugh
end_clip = 0;
start_clip = 0;
tc=3
figure(); hold on; title(strcat("M ",num2str(tc))); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]); 
scatter3(human_r(Mstart(tc)+start_clip:Mend(tc)-end_clip,1,4),human_r(Mstart(tc)+start_clip:Mend(tc)-end_clip,2,4),human_r(Mstart(tc)+start_clip:Mend(tc)-end_clip,3,4));



end