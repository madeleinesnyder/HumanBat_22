%% Trajectory extraction
x_mean = [ AllFlights(:,1)  AllFlights(:,2)  AllFlights(:,3)]'./1000;  

% Filter and interpolate
x_filt = x_mean; %medfilt1(x_mean,VideoFrameRate/2,[],2,'omitnan','truncate'); %filter after interpolating
x_intr = x_mean;% fillmissing(x_filt,'next',2,'EndValues','nearest');
x_spl = x_mean; %x_intr; %csaps(t, x_intr, 0.9, t);

%Frame rate
new_t = [1:length(AnalogSignals)];%t';
tracking_Fs = 120;VideoFrameRate;

%% Calculate velocity and flight segmentation
v = diff(x_spl,1,2)./[diff(new_t);diff(new_t);diff(new_t)]; v=[zeros(3,1) v];
v_abs = vecnorm(v,2,1);

nonflying = find(v_abs < 1);        toofast = find(v_abs > 30);
x_flying = x_spl;                   x_flying(:,[nonflying toofast]) = nan;
batspeed = v_abs;                   batspeed([nonflying toofast]) = nan;
bflying=~isnan(batspeed)';           %vector of 1s when the bat is flying

% For each sample, sum up the next 1s of data(flights are longer than 1s),Code from Nick D.
allsums = [];
for bf = 1 : size(bflying,1)-tracking_Fs
    allsums(bf) = sum(bflying(bf:bf+tracking_Fs));
end

% Detect flight starts and stops
[R,rLT,rUT,rLL,rUL] = risetime(allsums);    
[F,fLT,fUT,fLL,fUL] = falltime(allsums);           
if length(R) ~= length(F)
fLT(length(R)) = length(allsums);
fUT(length(R)) = length(allsums);
F(length(R)) = F(length(F));
end
flight_starts = round(rLT+tracking_Fs/2);
flight_ends = round(fLT+tracking_Fs/2); %... +Fs is a sistematic correction, useful
num_flights = size(R,2);
ref = ones(size(flight_starts));
avg_flight_time = mean((flight_ends-flight_starts)./tracking_Fs);

%% Plot results

% Plot 2D flight trajectories
plotFlightPathsAll = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
plot(x_mean(1,:),x_mean(2,:),'.');
hold on;        %rectangle('Position',[xL yB xR-xL yF-yB]);
%scatter([F3(1) F4(1)],[F3(2) F4(2)],'filled');  hold off;
xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
title(['Raw flights']);
xlabel('m'); ylabel('m');
hold off

subplot(1,2,2);
plot(x_spl(1,:),x_spl(2,:)); hold on; plot(x_mean(1,:),x_mean(2,:),'.','MarkerSize',1);
%rectangle('Position',[xL yB xR-xL yF-yB]);  hold off;
xlim([-3 3]); ylim([-3 3]); zlim([0 2.5]);
title(['Spline flights: ']);
xlabel('m'); ylabel('m');
hold off

% Plot session timeline
t=[1:length(AnalogSignals)];
plotFlightTimeline = figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(3,1,1);   plot(t,x_mean(1,:),'.');  hold on;
plot(new_t,x_spl(1,:),'--');   refline(0,F3(1));    hold off;
legend('cluster/mean','spl');   ylabel('x (m)');
ax2 = subplot(3,1,2);   plot(new_t,v_abs,'.');
hold on;    stem(new_t(flight_starts),ref);    stem(new_t(flight_ends),ref);  hold off;
ylabel('v (m/s)');
ax3 = subplot(3,1,3);  % plot(t,rew_signal);
ylabel('Rewards');
linkaxes([ax1,ax2,ax3],'x');    xlabel('Samples');