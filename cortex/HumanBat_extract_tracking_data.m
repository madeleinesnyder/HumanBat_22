function HumanBat_extract_tracking_data(work_dir)
c3dList = dir([work_dir filesep 'Generated_C3D_files' filesep '*-Bat_Cluster.c3d']);
% make processed directory
processed_dir = [work_dir filesep 'Generated_C3D_files' filesep 'processed' filesep];
if ~isdir(processed_dir)
    mkdir(processed_dir);
end


%run through list of c3d files in that directory, convert to mat, and save
%to processed directory
for i = 1:length(c3dList)
    disp(sprintf('Extracting %s', c3dList(i).name));
    fileName = extractBefore(c3dList(i).name,'-Bat_Cluster.c3d'); %Bat
    %fileName = extractBefore(c3dList(i).name,'-Bat_Cluster.c3d'); %Bat
    dateSesh = datestr(datetime(c3dList(i).date),'yymmdd');
    batName = extractBefore(fileName,['_' dateSesh]);
    %sessionNum = fileName(end);
    %copy_dir = [extractBefore(pwd,batName) 'processed' filesep batName filesep dateSesh filesep];
    %if ~isdir(copy_dir)
    %    mkdir(copy_dir);
    %end
    [Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]=HumanBat_readC3D_analog([work_dir filesep 'Generated_C3D_files' filesep c3dList(i).name]); %convert file
    %plot ttl impulses to check they are linear and not missing ttl
    event_ttls = AnalogSignals(:,2);
    [R,LT,UT,LL,UL] = risetime(event_ttls,VideoFrameRate);
    figure
    plot(LT)
    title(fileName)
    ttl_risetime = LT;
    %save new mat file in both original directory and copied directory for processing
    idx = Markers == 0;
    Markers(idx) = NaN;
    avgMarkerPos = squeeze(nanmean(Markers,2));

    %------------Local2GlobalTime-----------------
    local_ttl_usec = ttl_risetime*1e6;
    global_ttl_interval_usec = 3e6;
    local_sample_t0_usec = 0;%1e6/VideoFrameRate; % t0 starts at first frame?
    local_sample_ts_usec = linspace(local_sample_t0_usec, local_sample_t0_usec+(length(AnalogSignals)-1)*(1e6/VideoFrameRate), length(AnalogSignals));
    global_sample_ts_usec = local2GlobalTime(local_ttl_usec, local_sample_ts_usec, 'global_ttl_interval_us', 3e6);
    %---------------------------------------------
    
    % Save data
    save([processed_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate', 'avgMarkerPos','global_sample_ts_usec');
    %save([copy_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate');
end