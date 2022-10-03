function [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v6_ciholas_LAND(s,pts1,pts2,ids1,ids2,clrs1,clrs2,interval,cv1,cv2,ptitle1,ptitle2,PSTH_flag,fig_flag,Mx_flights,Kx_flights,ciholas_flight_struct_resort,clus)
%% Function for plotting spike rasters
%== INPUTS:
% s:            column vector with spike times
% pts:          column vector with event times for alignment
% ids:          column vector with the group identifier of each alignment event
% clrs:         nx3 matrix with the colors associated with each event
% interval:     pair of scalars [a,b] indicating the interval around the event time
% p_title:      text for the title
% PSTH_flag:    flag for plotting or not the PSTH
% fig_flag:     flag for plotting or not the raster (only calculates PSTH)

%== OUTPUTS:
% PSTH:         smoothed PSTH  
% cntrs_PSTH    time points corresponding to bin centers

%=== Make sure inputs have the correct size
if ~iscolumn(s),   s=s';        end     
if ~iscolumn(pts1), pts1=pts1';    end;      if ~iscolumn(pts2), pts2=pts2';    end
if ~iscolumn(ids1), ids1=ids1';    end;        if ~iscolumn(ids2), ids2=ids2';    end
if size(clrs1,2)~=3,clrs1=clrs1';  end;         if size(clrs2,2)~=3,clrs2=clrs2';  end

%=== Make sure that events are sorted and grouped
[ids1,sorted_indexes1] = sort(ids1);
pts1 = pts1(sorted_indexes1);
clrs1 = clrs1(sorted_indexes1,:);
[ids2,sorted_indexes2] = sort(ids2);
pts2 = pts2(sorted_indexes2);
clrs2 = clrs2(sorted_indexes2,:);

%=== Group definitions
g_clr1 = unique(clrs1,'rows','stable');
[group_id1,first_event1,~]= unique(ids1);
group_numerosity1 = diff([first_event1; length(pts1)+1]);
last_event1 = [first_event1(2:end)-1; length(pts1)];
g_clr2 = unique(clrs2,'rows','stable');
[group_id2,first_event2,~]= unique(ids2);
group_numerosity2 = diff([first_event2; length(pts2)+1]);
last_event2 = [first_event2(2:end)-1; length(pts2)];

%=== If only two groups, use default colors
group_id = group_id1;
if length(group_id)==2
   g_clr(1,:)=0.6*[1 1 1];
   g_clr(2,:)=[0.30,0.74,0.93];
end

%=== Parameters and Init
dt_PS = 0.1;
edges_PSTH = [interval(1):dt_PS:interval(2)];
PSTH_cs = edges_PSTH(1:end-1)+diff(edges_PSTH)/2;
PSTH_mn = zeros(numel(group_id),numel(PSTH_cs));
PSTH_up = zeros(numel(group_id),numel(PSTH_cs));
PSTH_dn = zeros(numel(group_id),numel(PSTH_cs));
A_factor1 = median(group_numerosity1)/4;        A_factor2 = median(group_numerosity2)/4; 
jj = 1; counter1 = 1;    coord_x1 = [];   coord_y1 = [];

%=== Calculate coordinates for raster plot and PSTH for all flights
list1=[]; list2=[];
[pts_all,pia] = sort([pts1;pts2]); ids_all = ([ids1;ids2]); 
for i=1:length(pts1)
    list1 = [list1,find(pia==i)];
end
for i=length(pts1)+1:length(pts2)+length(pts1)
    list2 = [list2,find(pia==i)];% (1:length(pts1)); list2 = pia(length(pts1)+1:end);
end
for i = group_id'
    chunk_starts = pts_all(ids_all==i);
    cum_spikes = []; 
    tmp_spikes = zeros(length(chunk_starts),numel(edges_PSTH)-1);
    for m = 1:length(chunk_starts)
        spikes = [];
        t_ref = chunk_starts(m);
        spikes = s(s>t_ref+interval(1) & s<t_ref+interval(2))-t_ref;
        tmp_spikes(m,:) = histcounts(spikes,edges_PSTH);
        cum_spikes = [cum_spikes; spikes];
        coord_x1 = cat(2,coord_x1,[spikes';spikes']);
        coord_y1 = cat(2,coord_y1,[ones(size(spikes'))-counter1;zeros(size(spikes'))-counter1]);
        counter1 = counter1+1;
    end
    PSTH_mn_(jj,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dt_PS*group_numerosity1(jj)),'gaussian',5);
    PSTH_mn_(jj,:) = smoothdata(mean(tmp_spikes,1)/dt_PS,'gaussian',5);
    PSTH_up(jj,:) = smoothdata(prctile(tmp_spikes,95,1)/dt_PS,'gaussian',5);
    PSTH_dn(jj,:) = smoothdata(prctile(tmp_spikes, 5,1)/dt_PS,'gaussian',5);
    jj = jj +1;
end

% For the lines at the bottom of the plot, make seperate PSTH things
% For Madeleine flights
clear PSTH_up PSTH_dn PSTH_mn1
jj_1 = 1; counterfake1 = 1;    coord_fakex1 = [];   coord_fakey1 = [];
chunk_starts = pts1;
cum_spikes = []; 
tmp_spikes = zeros(length(chunk_starts),numel(edges_PSTH)-1);
for m = 1:length(chunk_starts)
    spikes = [];
    t_ref = chunk_starts(m);
    spikes = s(s>t_ref+interval(1) & s<t_ref+interval(2))-t_ref;
    tmp_spikes(m,:) = histcounts(spikes,edges_PSTH);
    cum_spikes = [cum_spikes; spikes];
    coord_fakex1 = cat(2,coord_fakex1,[spikes';spikes']);
    coord_fakey1 = cat(2,coord_fakey1,[ones(size(spikes'))-counterfake1;zeros(size(spikes'))-counterfake1]);
    counterfake1 = counterfake1+1;
end
PSTH_mn1(jj_1,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dt_PS*group_numerosity1(jj_1)),'gaussian',5);
PSTH_mn1(jj_1,:) = smoothdata(mean(tmp_spikes,1)/dt_PS,'gaussian',5);
PSTH_up(jj_1,:) = smoothdata(prctile(tmp_spikes,95,1)/dt_PS,'gaussian',5);
PSTH_dn(jj_1,:) = smoothdata(prctile(tmp_spikes, 5,1)/dt_PS,'gaussian',5);
jj_1 = jj_1 +1;

% For the kevin flights
jj_2 = 1; counterfake2 = 1;    coord_fakex2 = [];   coord_fakey2 = [];
chunk_starts = pts2;
cum_spikes = []; 
tmp_spikes = zeros(length(chunk_starts),numel(edges_PSTH)-1);
for m = 1:length(chunk_starts)
    spikes = [];
    t_ref = chunk_starts(m);
    spikes = s(s>t_ref+interval(1) & s<t_ref+interval(2))-t_ref;
    tmp_spikes(m,:) = histcounts(spikes,edges_PSTH);
    cum_spikes = [cum_spikes; spikes];
    coord_fakex2 = cat(2,coord_fakex2,[spikes';spikes']);
    coord_fakey2 = cat(2,coord_fakey2,[ones(size(spikes'))-counterfake2;zeros(size(spikes'))-counterfake2]);
    counterfake2 = counterfake2+1;
end
PSTH_mn2(jj_2,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dt_PS*group_numerosity1(jj_2)),'gaussian',5);
PSTH_mn2(jj_2,:) = smoothdata(mean(tmp_spikes,1)/dt_PS,'gaussian',5);
PSTH_up(jj_2,:) = smoothdata(prctile(tmp_spikes,95,1)/dt_PS,'gaussian',5);
PSTH_dn(jj_2,:) = smoothdata(prctile(tmp_spikes, 5,1)/dt_PS,'gaussian',5);
jj_2 = jj_2 +1;

last_event_all = [first_event1(2:end)-1; length(pts_all)];

%=== Plot raster
cluster_cols1 = jet(length(unique(cv1)));
cluster_cols = cluster_cols1;

if fig_flag
    figure();
    subplot(1,2,1); hold on; 
    for jj=1:length(coord_y1)
        if ismember(-coord_y1(2,jj),list1)
            plot(coord_x1(:,jj),coord_y1(:,jj),'b-','LineWidth',1);    %hold on; 
        elseif ismember(-coord_y1(2,jj),list2)
            plot(coord_x1(:,jj),coord_y1(:,jj),'r-','LineWidth',1);    %hold on; 
        end
    end
    for j = 1:numel(group_id)
        line([interval(1),interval(2)], [-last_event_all(j),-last_event_all(j)],'Color', g_clr1(j,:),'LineStyle','-','LineWidth',2);
        if PSTH_flag
            plot(PSTH_cs,(A_factor1/max(PSTH_mn1(:)))*PSTH_mn1(j,:)-last_event_all(j),'Color', [0.48 0.72 0.98],'LineWidth',3);
            plot(PSTH_cs,(A_factor2/max(PSTH_mn2(:)))*PSTH_mn2(j,:)-last_event_all(j),'Color', [1 0.5 0.8],'LineWidth',3);
        end
    end
    line([0,0], [0,-counter1+1],'Color', [1 0 0],'LineStyle','--');    %hold off;
    ylim([-counter1+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
    title(ptitle1);
    h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
    hold off;

    subplot(1,2,2); hold on; title(strcat("All Cluster ",num2str(clus)," Flights")); xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    K_cols = [1 0 0 1; 1 0 0 0.8; 1 0 0 0.6; 1 0 0 0.4; 1 0 0 0.2];
    M_cols = [0 0 1 1; 0 0 1 0.8; 0 0 1 0.6; 0 0 1 0.4; 0 0 1 0.2];
    K_cols = [1 0 0; 1 0.2 0.2; 1 0.3 0.3; 1 0.4 0.4; 1 0.5 0.5; 1 0.55 0.55; 1 0.6 0.6; 1 0.65 0.65; 1 0.7 0.7; 1 0.75 0.75; 1 0.8 0.8;];
    M_cols = [0 0 1 ; 0.2 0.2 1; 0.3 0.3 1; 0.4 0.4 1; 0.5 0.5 1; 1 0.55 0.55; 1 0.6 0.6; 1 0.65 0.65; 1 0.7 0.7; 1 0.75 0.75; 1 0.8 0.8;];
    for mm = 1:length(Mx_flights)
            plot3(ciholas_flight_struct_resort{Mx_flights(mm)}.pos(:,1),ciholas_flight_struct_resort{Mx_flights(mm)}.pos(:,2),ciholas_flight_struct_resort{Mx_flights(mm)}.pos(:,3),'Color',M_cols(mm,:)); 
            scatter3(ciholas_flight_struct_resort{Mx_flights(mm)}.pos(1,1),ciholas_flight_struct_resort{Mx_flights(mm)}.pos(1,2),ciholas_flight_struct_resort{Mx_flights(mm)}.pos(1,3),50,'k','filled');
    end
    for mm = 1:length(Kx_flights)
            plot3(ciholas_flight_struct_resort{Kx_flights(mm)}.pos(:,1),ciholas_flight_struct_resort{Kx_flights(mm)}.pos(:,2),ciholas_flight_struct_resort{Kx_flights(mm)}.pos(:,3),'Color',K_cols(mm,:)); 
            scatter3(ciholas_flight_struct_resort{Kx_flights(mm)}.pos(1,1),ciholas_flight_struct_resort{Kx_flights(mm)}.pos(1,2),ciholas_flight_struct_resort{Kx_flights(mm)}.pos(1,3),50,'k','filled');
    end
    hold off;
end

end