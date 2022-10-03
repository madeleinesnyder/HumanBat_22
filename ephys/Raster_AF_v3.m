function [PSTH_mn,PSTH_cs,PSTH_up,PSTH_dn] = Raster_AF_v3(s,pts,ids,clrs,interval,cv,p_title,PSTH_flag,fig_flag)
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
if ~iscolumn(pts), pts=pts';    end
if ~iscolumn(ids), ids=ids';    end
if size(clrs,2)~=3,clrs=clrs';  end

%=== Make sure that events are sorted and grouped
[ids,sorted_indexes] = sort(ids);
pts = pts(sorted_indexes);
clrs = clrs(sorted_indexes,:);

%=== Group definitions
g_clr = unique(clrs,'rows','stable');
[group_id,first_event,~]= unique(ids);
group_numerosity = diff([first_event; length(pts)+1]);
last_event = [first_event(2:end)-1; length(pts)];

%=== If only two groups, use default colors
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
A_factor = median(group_numerosity)/4; 
jj = 1; counter = 1;    coord_x = [];   coord_y = [];

%=== Calculate coordinates for raster plot and PSTH
for i = group_id'
    chunk_starts = pts(ids==i);
    cum_spikes = []; 
    tmp_spikes = zeros(length(chunk_starts),numel(edges_PSTH)-1);
    for m = 1:length(chunk_starts)
        spikes = [];
        t_ref = chunk_starts(m);
        spikes = s(s>t_ref+interval(1) & s<t_ref+interval(2))-t_ref;
        tmp_spikes(m,:) = histcounts(spikes,edges_PSTH);
        cum_spikes = [cum_spikes; spikes];
        coord_x = cat(2,coord_x,[spikes';spikes']);
        coord_y = cat(2,coord_y,[ones(size(spikes'))-counter;zeros(size(spikes'))-counter]);
        counter = counter+1;
    end
    PSTH_mn(jj,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dt_PS*group_numerosity(jj)),'gaussian',5);
    PSTH_mn(jj,:) = smoothdata(mean(tmp_spikes,1)/dt_PS,'gaussian',5);
    PSTH_up(jj,:) = smoothdata(prctile(tmp_spikes,95,1)/dt_PS,'gaussian',5);
    PSTH_dn(jj,:) = smoothdata(prctile(tmp_spikes, 5,1)/dt_PS,'gaussian',5);
    jj = jj +1;
end

%=== Plot raster
cluster_cols = jet(length(unique(cv)));
if ~isempty(cv)
    if fig_flag
        h = figure(); hold on;
        plot(coord_x,coord_y,'k-','LineWidth',1);    %hold on; 
        for j = 1:numel(group_id)
            line([interval(1),interval(2)], [-last_event(j),-last_event(j)],'Color', g_clr(j,:),'LineStyle','-','LineWidth',2);
            if PSTH_flag
                plot(PSTH_cs,(A_factor/max(PSTH_mn(:)))*PSTH_mn(j,:)-last_event(j),'Color', [g_clr(j,:) 0.7],'LineWidth',2);
            end
        end
        line([0,0], [0,-counter+1],'Color', 0.5*[1 1 1],'LineStyle','--');    %hold off;
        %figure(); ylim([-counter+1 0]); xlim(interval);
        u_cv = unique(cv);
        for j=1:length(unique(cv))
            cluster_limits_ = find(cv==u_cv(j));
            line([interval(1) interval(2)], [-cluster_limits_(end) -cluster_limits_(end)],'Color',cluster_cols(j,:),'LineStyle','--');
        end
        ylim([-counter+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
        title(p_title);
        %h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        %h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        hold off;
    end
else
    if fig_flag
        figure(); hold on;
        plot(coord_x,coord_y,'k-','LineWidth',1);    %hold on; 
        for j = 1:numel(group_id)
            line([interval(1),interval(2)], [-last_event(j),-last_event(j)],'Color', g_clr(j,:),'LineStyle','-','LineWidth',2);
            if PSTH_flag
                plot(PSTH_cs,(A_factor/max(PSTH_mn(:)))*PSTH_mn(j,:)-last_event(j),'Color', [g_clr(j,:) 0.7],'LineWidth',2);
            end
        end
        line([0,0], [0,-counter+1],'Color', 0.5*[1 1 1],'LineStyle','--');    %hold off;
        ylim([-counter+1 0]); xlim(interval); ylabel('Trial #'); xlabel('Time (s)');
        title(p_title);
        h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual'); %h.YTickLabelMode = 'manual';  % Correct trial sign
        hold off;
    end
end
end