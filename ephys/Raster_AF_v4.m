function [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v4(s,pts1,pts2,ids1,ids2,clrs1,clrs2,interval,cv1,cv2,ptitle1,ptitle2,PSTH_flag,fig_flag)
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

%=== Calculate coordinates for raster plot and PSTH #1
for i = group_id'
    chunk_starts = pts1(ids1==i);
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
    PSTH_mn1(jj,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dt_PS*group_numerosity1(jj)),'gaussian',5);
    PSTH_mn1(jj,:) = smoothdata(mean(tmp_spikes,1)/dt_PS,'gaussian',5);
    PSTH_up(jj,:) = smoothdata(prctile(tmp_spikes,95,1)/dt_PS,'gaussian',5);
    PSTH_dn(jj,:) = smoothdata(prctile(tmp_spikes, 5,1)/dt_PS,'gaussian',5);
    jj = jj +1;
end

%=== Calculate coordinates for raster plot and PSTH #2
jj=1; counter2 = 1;    coord_x2 = [];   coord_y2 = [];

for i = group_id'
    chunk_starts = pts2(ids2==i);
    cum_spikes = []; 
    tmp_spikes = zeros(length(chunk_starts),numel(edges_PSTH)-1);
    for m = 1:length(chunk_starts)
        spikes = [];
        t_ref = chunk_starts(m);
        spikes = s(s>t_ref+interval(1) & s<t_ref+interval(2))-t_ref;
        tmp_spikes(m,:) = histcounts(spikes,edges_PSTH);
        cum_spikes = [cum_spikes; spikes];
        coord_x2 = cat(2,coord_x2,[spikes';spikes']);
        coord_y2 = cat(2,coord_y2,[ones(size(spikes'))-counter2;zeros(size(spikes'))-counter2]);
        counter2 = counter2+1;
    end
    PSTH_mn2(jj,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dt_PS*group_numerosity2(jj)),'gaussian',5);
    PSTH_mn2(jj,:) = smoothdata(mean(tmp_spikes,1)/dt_PS,'gaussian',5);
    PSTH_up(jj,:) = smoothdata(prctile(tmp_spikes,95,1)/dt_PS,'gaussian',5);
    PSTH_dn(jj,:) = smoothdata(prctile(tmp_spikes, 5,1)/dt_PS,'gaussian',5);
    jj = jj +1;
end

%=== Plot raster
cluster_cols1 = jet(length(unique(cv1)));
cluster_cols2 = jet(length(unique(cv2)));
if size(cluster_cols2,1) > size(cluster_cols1,1)
    cluster_cols = cluster_cols2;
elseif size(cluster_cols1,1) > size(cluster_cols2,1)
    cluster_cols = cluster_cols1;
else 
    cluster_cols = cluster_cols2;
end
if ~isempty(cv1)
    if fig_flag

        h = figure(); subplot(1,2,1); hold on;
        plot(coord_x1,coord_y1,'k-','LineWidth',1);    %hold on; 
        for j = 1:numel(group_id)
            p1=line([interval(1),interval(2)], [-last_event1(j),-last_event1(j)],'Color', g_clr1(j,:),'LineStyle','-','LineWidth',2);
            if PSTH_flag
                plot(PSTH_cs,(A_factor1/max(PSTH_mn1(:)))*PSTH_mn1(j,:)-last_event1(j),'Color', [g_clr1(j,:) 0.7],'LineWidth',2);
            end
        end
        cv1_legend = string(unique(cv1)); 
        %legend("","","","","","","");
        line([0,0], [0,-counter1+1],'Color', [1 0 0],'LineStyle','--');    %hold off;
        %figure(); ylim([-counter+1 0]); xlim(interval);
        u_cv1 = unique(cv1);
        for j=1:length(unique(cv1))
            cluster_limits_ = find(cv1==u_cv1(j));
            p2(j)=line([interval(1) interval(2)], [-cluster_limits_(end) -cluster_limits_(end)],'Color',cluster_cols(j,:),'LineStyle','--');
        end
        ylim([-counter1+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
        title(ptitle1);
        legend(p2,cv1_legend);
        %h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        %h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        hold off;

        subplot(1,2,2); hold on;
        plot(coord_x2,coord_y2,'k-','LineWidth',1);    %hold on; 
        for j = 1:numel(group_id)
            p1 = line([interval(1),interval(2)], [-last_event2(j),-last_event2(j)],'Color', g_clr2(j,:),'LineStyle','-','LineWidth',2);
            if PSTH_flag
                plot(PSTH_cs,(A_factor2/max(PSTH_mn2(:)))*PSTH_mn2(j,:)-last_event2(j),'Color', [g_clr2(j,:) 0.7],'LineWidth',2);
            end
        end
        line([0,0], [0,-counter2+1],'Color', [1 0 0],'LineStyle','--');    %hold off;
        %figure(); ylim([-counter+1 0]); xlim(interval);
        u_cv2 = unique(cv2);
        cv2_legend = string(u_cv2);
        for j=1:length(unique(cv2))
            cluster_limits_ = find(cv2==u_cv2(j));
            p2(j) = line([interval(1) interval(2)], [-cluster_limits_(end) -cluster_limits_(end)],'Color',cluster_cols(j,:),'LineStyle','--');
        end
        ylim([-counter2+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
        title(ptitle2);
        legend(p2,cv2_legend);
        %h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        %h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        hold off;
    end
else
    if fig_flag
        figure();
        subplot(1,2,1); hold on; plot(coord_x1,coord_y1,'k-','LineWidth',1);    %hold on; 
        for j = 1:numel(group_id)
            line([interval(1),interval(2)], [-last_event1(j),-last_event1(j)],'Color', g_clr1(j,:),'LineStyle','-','LineWidth',2);
            if PSTH_flag
                plot(PSTH_cs,(A_factor1/max(PSTH_mn1(:)))*PSTH_mn1(j,:)-last_event1(j),'Color', [g_clr1(j,:) 0.7],'LineWidth',2);
            end
        end
        line([0,0], [0,-counter1+1],'Color', [1 0 0],'LineStyle','--');    %hold off;
        ylim([-counter1+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
        title(ptitle1);
        h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        hold off;

        subplot(1,2,2); hold on; hold on; plot(coord_x2,coord_y2,'k-','LineWidth',1);    %hold on; 
        for j = 1:numel(group_id)
            line([interval(1),interval(2)], [-last_event2(j),-last_event2(j)],'Color', g_clr2(j,:),'LineStyle','-','LineWidth',2);
            if PSTH_flag
                plot(PSTH_cs,(A_factor2/max(PSTH_mn2(:)))*PSTH_mn2(j,:)-last_event2(j),'Color', [g_clr2(j,:) 0.7],'LineWidth',2);
            end
        end
        line([0,0], [0,-counter2+1],'Color', [1 0 0],'LineStyle','--');    %hold off;
        ylim([-counter2+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
        title(ptitle2);
        h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        hold off;
    end
end
end