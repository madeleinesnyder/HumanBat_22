function [outputArg1,outputArg2] = SplitCdpChunks(fname,varargin)
%SPLITCDPCHUNKS Split CDP data into chunks based on TTLs (e.g. rest,fly,rest)
%   Splits CDP into chunks denoted by chunks of TTLs
%   Output is to be then processed by ExtractCdp_AF_v0.m
%-----------------------------------------------------------------
% Data are recorded through a python script and saved in a ASCII file
% Data columns (comma delimited) are defines as follows:

% Column 1 = pos(1), sync(2) or acc(0) type index
% Column 2 = device Serial Number
% Column 3 = Network time (each tic is ~15.65 ps)
% Column 4,5,6 = x,y,z (position in mm or acceleration-sync)
% Column 7 = signal quality for pos or scale for acc
% Column 8 = number of receiving anchors for a tag or '0' for acc/sync

% The sync signal is driven by a TTL:
% 50 ms duration happening every 21, 13, 8, 5, 4s
%-----------------------------------------------------------------
%   *********USAGE EXAMPLES*****************
% ExtractCdp_AF_v0('_cdp_fly_1','Exclude',17106999);
% ExtractCdp_AF_v0('_cdp_fly_2','Include',[17106917 17106934 17107055 17106969 17106904 17106951]);
% ExtractCdp_AF_v0('_cdp_fly_1','Include',[17106917 17106934 17107055 17106969 17106904 17106951 17107468]);

%Load data, keep non-zero entries, sort according to ascending network times and delete duplicates
file_name = fname
disp(file_name)

parts = strsplit(file_name, filesep);
parent_path = strjoin(parts(1:end-1), filesep);
fn = parts(end);
fn = fn{1}(1:end-4);
disp([parent_path, '/extracted_', fn, '.mat'])

RTLS_data = load(file_name);
RTLS_data = RTLS_data(RTLS_data(:,2)~=0,:);
RTLS_data = sortrows(RTLS_data,3);
RTLS_data = unique(RTLS_data,'rows','stable');

%User inputs overrides
tags_SN = [17106917; 17106934; 17107055; 17106969; 17106904; 17107000; 17106999; 17106951;  17107429;   17107383];
if nargin > 1
    eSN = [];   iSN = [];
    nparams=length(varargin);
    for i=1:2:nparams
        switch (varargin{i})
            case 'Exclude'
                eSN=varargin{i+1};
                tags_SN(any(tags_SN == eSN,2)) = [];
            case 'Include'
                iSN = varargin{i+1};
                tags_SN = [];
                
                if iscell(iSN')
                    tags_SN = cell2mat(iSN)';
                else
                    tags_SN = iSN';
                end
        end
    end
end
disp(tags_SN)
%Parameters
use_sync = 1;
rec_duration = 10800;                                                                           %approx rec duration (s), 10800 covers 3 hours
ntw_time = 15.65e-12;                                                                           %network tic interval (s)
%tags_SN = [17106917; 17106934; 17107055; 17106969; 17107100; 17106999];                        %tag serial numbers
%tags_SN = [17106917; 17106934; 17107055; 17106969; 17106904; 17106951];                        %tag serial numbers
sync_SN = 17040920; %17106963;                                                                  %sync tag serial number
n_tags = length(tags_SN);
TTL_time_diff = [21; 13; 8; 5; 4];                                                              %TTL delays in s
TTL_abs_times = [0; cumsum(repmat(TTL_time_diff,round(rec_duration*2/sum(TTL_time_diff)),1))];
CDPmtdata.Fs = 100;                                                                             %Acquisition frequency CDP (Hz)
CDPmtdata.tag_SN = tags_SN;
CDPmtdata.sync_SN = sync_SN;

%Extract data and log CDP metadata, start from sync
sync_data = RTLS_data(RTLS_data(:,1)==2 & RTLS_data(:,2)==sync_SN,2:end);
CDPmtdata.sync = ~isempty(sync_data);
CDPmtdata.sync_duration = (sync_data(end,2)-sync_data(1,2))*ntw_time/60;
CDPmtdata.sync_samples = length(sync_data);

%Tag Data
j = 0;
for i = 1:n_tags
    disp(tags_SN(i))
    if ~isempty(find(RTLS_data(:,2)==tags_SN(i)))
        j = j+1;
        tag_data{j} = RTLS_data(RTLS_data(:,1)==1 & RTLS_data(:,2)==tags_SN(i),2:end);
        tag_data{1,j}(:,[3:5]) = tag_data{1,j}(:,[3:5])/1000;
        CDPmtdata.tag_duration(j) = (tag_data{1,j}(end,2)-tag_data{1,j}(1,2))*ntw_time/60;
        CDPmtdata.tag_samples(j) = length(tag_data{1,j});
        
        if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(i))))
            tag_ac_data{j} = RTLS_data(RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(i),2:end);
            %Acceleration values need to be converted
            indexes = [false(size(tag_ac_data{1,j}(:,1:2))), tag_ac_data{1,j}(:,[3:5])>double(intmax('int32'))];
            tag_ac_data{1,j}(indexes) = tag_ac_data{1,j}(indexes)-double(intmax('uint32'));
            tag_ac_data{1,j}(:,[3:5]) = tag_ac_data{1,j}(:,[3:5])*2/double(intmax('int32'));
            CDPmtdata.tag_ac_duration(j) = (tag_ac_data{1,j}(end,2)-tag_ac_data{1,j}(1,2))*ntw_time/60;
            CDPmtdata.tag_ac_samples(j) = length(tag_ac_data{1,j});
        end
    end
end

CDPmtdata.tags = n_tags;
disp(CDPmtdata);

[~, TTL_network_times] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time);

end

