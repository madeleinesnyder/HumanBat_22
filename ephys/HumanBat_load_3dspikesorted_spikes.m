%% Load Sorted Spikes (.ntt) and create cell structures with all the relevant information

%=== Load data and set parameters
TT_unit_names = dir(fullfile(cd, '*_SS_*.ntt'));
batdate = TT_unit_names(1).name(regexp(TT_unit_names(1).name,'_2022')+[3:8]);
load('TTL_timestamps');

n_units = length(TT_unit_names);
TT_unit = struct([]);
uV_factor = 0.0305; %or 3.3uV?
n_samples = 32;
Fs = 31.25e3;
t = [1:n_samples]/Fs; 
save_data = 1;
plot_fig = 1;
savefigures = 1;  
fig_count = 1;

if save_data
    analysis_directory=fullfile(pwd,['SingleUnits_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end  
end

%% Process each unit

for i = 1:n_units
    
    %=== Extract Single Unit Data
    [TT_unit(i).Timestamps, ~, TT_unit(i).CellNumbers, TT_unit(i).Features, TT_unit(i).Samples, TT_unit(i).Header] = Nlx2MatSpike(TT_unit_names(i).name, [1 1 1 1 1], 1, 1, [] );
    TT_unit(i).Samples = TT_unit(i).Samples*uV_factor;
    TT_unit(i).TT = str2double(TT_unit_names(i).name(regexp(TT_unit_names(1).name,'_TT')+3));
    TT_unit(i).NSpikes = length(TT_unit(i).Timestamps);
    TT_unit(i).Refr_viol = nnz(diff(TT_unit(i).Timestamps)<1e3)/TT_unit(i).NSpikes;
    TT_unit(i).Avg_wform = squeeze(mean(TT_unit(i).Samples(:,:,:),3));
    TT_unit(i).Frng_freq = 1e6/mean(diff(TT_unit(i).Timestamps));
    TT_unit(i).Isdist = 1;
    TT_unit(i).Lratio = 1;
    TT_unit(i).Fs = Fs;
    TT_unit(i).amp = max(abs(max(TT_unit(i).Avg_wform)-min(TT_unit(i).Avg_wform)));
    TT_unit(i).Incl_flag = 1;
    
    %=== Show features of the unit
    if plot_fig
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        %=== Plot a subset of spike waveforms
        subplot(5,4,[1,2,5,6]);
        if TT_unit(i).NSpikes>=1e3, subset = datasample(1:TT_unit(i).NSpikes,1e3,'Replace',false);
        else, subset = 1:TT_unit(i).NSpikes;end
        Subset_samples = TT_unit(i).Samples(:,:,subset);
        plot(reshape(cat(1,Subset_samples,nan(1,4,size(Subset_samples,3))),33*4,[]),'k'); hold on;
        plot(reshape([TT_unit(i).Avg_wform; nan(1,4)],[],1),'r','LineWidth',3);  hold off;
        ylim([min(Subset_samples,[],'all') max(Subset_samples,[],'all')]);
        xlabel('samples');    ylabel('uV');
        sgtitle(['TT' num2str(TT_unit(i).TT) ' Unit' num2str(TT_unit(i).CellNumbers(1)) ': ' num2str(TT_unit(i).NSpikes) ' spikes']);
        
        %=== Plot ISI histogram
        subplot(5,4,[3,4,7,8]);
        edges_us = 10.^(0:0.05:9);
        histogram(diff(TT_unit(i).Timestamps),edges_us,'edgecolor','none');
        set(gca,'XScale','log');
        y1=get(gca,'ylim');  hold on; plot([1e3 1e3],y1);   hold off;
        xlabel('us');    ylabel('ISI count');
        title(['Average Firing Frequency: ' num2str(TT_unit(i).Frng_freq,3) ' Hz']);
        
        %=== Plot amplitude (on larger channel) throughout the session
        ax(1) = subplot(5,4,[9:10,13:14]);
        [maxVal,maxCh] = max(max(TT_unit(i).Avg_wform));
        [~,maxSmp] = max(TT_unit(i).Avg_wform(:,maxCh));
        plot(TT_unit(i).Timestamps, squeeze(TT_unit(i).Samples(maxSmp,maxCh,:)));
        hold on; refline(0,maxVal); y1=get(gca,'ylim');
        plot([TTL_timestamps(1) TTL_timestamps(1)],y1); hold off;
        % plot([TTL_timestamps.pre(1) TTL_timestamps.pre(1)],y1,[TTL_timestamps.fly(1) TTL_timestamps.fly(1)],y1,[TTL_timestamps.pst(1) TTL_timestamps.pst(1)],y1); hold off;
        xlabel('Timestamp');    ylabel('Spike Amplitude');
        
        %=== Plot spike counts (~firing freq) in 1s ms windows
        ax(2) = subplot(5,4,[17:18]);
        edges_fr = [TT_unit(i).Timestamps(1)-1e6:1e6:TT_unit(i).Timestamps(end)+1e6];
        histogram(TT_unit(i).Timestamps,edges_fr,'edgecolor','none','FaceColor','black');
        xlabel('Timestamp');    ylabel('Spike count');
        linkaxes(ax,'x');
        
        %=== Plot waveform on thirds
        dl = (TT_unit(i).Timestamps(end)-TT_unit(i).Timestamps(1))/3;
        first_3rd = TT_unit(i).Samples(:,:,find(TT_unit(i).Timestamps>TT_unit(i).Timestamps(1)+0*dl & TT_unit(i).Timestamps<TT_unit(i).Timestamps(1)+1*dl));
        secnd_3rd = TT_unit(i).Samples(:,:,find(TT_unit(i).Timestamps>TT_unit(i).Timestamps(1)+1*dl & TT_unit(i).Timestamps<TT_unit(i).Timestamps(1)+2*dl));
        third_3rd = TT_unit(i).Samples(:,:,find(TT_unit(i).Timestamps>TT_unit(i).Timestamps(1)+2*dl & TT_unit(i).Timestamps<TT_unit(i).Timestamps(1)+3*dl));
        subplot(5,4,[11,12,15,16,19,20]);
        plot(reshape([mean(first_3rd,3); ones(1,4)],[],1),'r');  hold on;
        ciplot(reshape([mean(first_3rd,3); ones(1,4)],[],1)-reshape([std(first_3rd,[],3); ones(1,4)],[],1),...
            reshape([mean(first_3rd,3); ones(1,4)],[],1)+reshape([std(first_3rd,[],3); ones(1,4)],[],1),...
            1:size(reshape([mean(first_3rd,3); NaN(1,4)],[],1)),'r');   alpha(0.3);
        plot(reshape([mean(secnd_3rd,3); ones(1,4)],[],1),'g');  
        ciplot(reshape([mean(secnd_3rd,3); ones(1,4)],[],1)-reshape([std(secnd_3rd,[],3); ones(1,4)],[],1),...
            reshape([mean(secnd_3rd,3); ones(1,4)],[],1)+reshape([std(secnd_3rd,[],3); ones(1,4)],[],1),...
            1:size(reshape([mean(secnd_3rd,3); NaN(1,4)],[],1)),'g');   alpha(0.3);
        plot(reshape([mean(third_3rd,3); ones(1,4)],[],1),'b');  
        ciplot(reshape([mean(third_3rd,3); ones(1,4)],[],1)-reshape([std(third_3rd,[],3); ones(1,4)],[],1),...
            reshape([mean(third_3rd,3); ones(1,4)],[],1)+reshape([std(third_3rd,[],3); ones(1,4)],[],1),...
            1:size(reshape([mean(third_3rd,3); NaN(1,4)],[],1)),'b');   alpha(0.3); hold off;
        legend('1st Third','','2nd Third','','3rd Third');  ylabel('uV');
        
        %=== Look at cluster changing across session
        %... Work in progress...
    end
    
    answer = questdlg('How would you score the unit?','Single Unit Inspection','Good','Uncertain','Bad','Good');
    switch answer
        case 'Good'
            TT_unit(i).Incl_flag = 1;
        case 'Uncertain'
            TT_unit(i).Incl_flag = NaN;
        case 'Bad'
            TT_unit(i).Incl_flag = 0;
    end
    
    fig_count = saveFig(analysis_directory,batdate,fig_count,savefigures);
    
    
end

%===Save data
if save_data
    save([analysis_directory,'/SingleUnits_', batdate, '.mat'],'TT_unit');
end