function PSTH = Raster_TimeWrap_AF_v1(s,pts_1,pts_2,ids,clrs,f_prc,p_title,PSTH_flag)
%% Function for plotting spike rasters by warping time of the event
%== INPUTS:
% s:            column vector with spike times
% pts:          column vector with event start times for alignment
% pts:          column vector with event stop  times for alignment
% ids:          column vector with the group identifier of each alignment event
% clrs:         nx3 matrix with the colors associated with each event
% f_prc:        percent of the event duration to be inspected before and after
% p_title:      text for the title
% PSTH_flag:    flag for plotting or not the PSTH
% fig_flag:     flag for plotting or not the raster (only calculates PSTH)

%== OUTPUTS:
% PSTH:         smoothed PSTH  
% cntrs_PSTH    time points corresponding to bin centers

%=== Make sure inputs have the correct size
if ~iscolumn(s),   s=s';              end
if ~iscolumn(pts_1), pts_1=pts_1';    end
if ~iscolumn(pts_2), pts_2=pts_2';    end
if ~iscolumn(ids), ids=ids';          end
if size(clrs,2)~=3,clrs=clrs';        end

%=== Make sure that events are sorted and grouped
[ids,sorted_indexes] = sort(ids);
pts_1 = pts_1(sorted_indexes);
pts_2 = pts_2(sorted_indexes);
clrs = clrs(sorted_indexes,:);

%=== Group definitions
g_clr = unique(clrs,'rows','stable');
[group_id,first_event,~]= unique(ids);
group_numerosity = diff([first_event; length(pts_1)+1]);
last_event = [first_event(2:end)-1; length(pts_1)];

%=== If only two groups, use default colors
if length(group_id)==2
   g_clr(1,:)=0.6*[1 1 1];
   g_clr(2,:)=[0.30,0.74,0.93];
end

% Params and init
dp_PS = 0.05;            %Percentage
edges_PSTH = [-f_prc:dp_PS:1+f_prc];
cntrs_PSTH = edges_PSTH(1:end-1)+diff(edges_PSTH)/2;
A_factor = median(group_numerosity)/4; 
PSTH = zeros(numel(group_id),numel(cntrs_PSTH));
jj = 1; counter = 1;    coord_x = [];   coord_y = [];

%=== Calculate coordinates for raster plot and PSTH
for i = group_id'
    chunk_strts = pts_1(ids==i);
    chunk_stops = pts_2(ids==i);
    cum_spikes = [];  
    for m = 1:length(chunk_strts)
        spikes = [];
        t_ref1 = chunk_strts(m);
        t_ref2 = chunk_stops(m);
        dur = t_ref2-t_ref1;
        spikes = s(s>t_ref1-dur*f_prc & s<t_ref2+dur*f_prc)-t_ref1; 
        spikes_int_real = interp1(linspace(-dur*f_prc,dur+dur*f_prc,1000*(1+2*f_prc)),linspace(-f_prc,1+f_prc,1000*(1+2*f_prc)),spikes);
        cum_spikes = [cum_spikes; spikes_int_real];
        coord_x = cat(2,coord_x,[spikes_int_real';spikes_int_real']);
        coord_y = cat(2,coord_y,[ones(size(spikes_int_real'))-counter;zeros(size(spikes_int_real'))-counter]);
        counter = counter+1;
    end
    PSTH(jj,:) = smoothdata(histcounts(cum_spikes,edges_PSTH)/(dp_PS*group_numerosity(jj)),'gaussian',5);
    jj = jj +1;
end

%=== Plot raster
plot(coord_x,coord_y,'k-','LineWidth',1);    hold on;
for j = 1:numel(group_id)
    line([-f_prc,1+f_prc], [-last_event(j),-last_event(j)],'Color', g_clr(j,:),'LineStyle','-','LineWidth',2);
    if PSTH_flag
        plot(edges_PSTH(1:end-1),(A_factor/max(PSTH(:)))*PSTH(j,:)-last_event(j),'Color', [g_clr(j,:) 0.7],'LineWidth',2);
    end
end
line([0,0], [0,-counter+1],'Color', 0.5*[1 1 1],'LineStyle','--');
line([1,1], [0,-counter+1],'Color', 0.7*[1 1 1],'LineStyle',':');
ylim([-counter+1 0]); xlim([-f_prc,1+f_prc]); ylabel('Trial #'); xlabel('Flight %');
title(p_title);
h = gca;    h.YTickLabel = num2str(abs(h.YTick)');  yticklabels('manual'); yticks('manual');% Correct trial sign
hold off;

end

