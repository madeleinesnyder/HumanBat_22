function [alignedCiholasCortex] = HumanBat_alignCiholasCortex(cortex, ciholas, ciholas2cortex)
%alignCiholasCortex Align Ciholas and Cortex data: 
%       1) Same position scales                                     
%       2) Same sampling frequency (cortex frequency, 120hz)
%       3) Same starting sample (sample at first TTL)
%   Parameters
%   ----------
%   cortex
%       extracted cortex data (contains Markers, AnalogSignals, etc...)
%   ciholas
%       extracted ciholas data from ExtractCdp_AF_v0.m
%   ciholas2cortex
%       Scaling factors to rescale ciholas to cortex units. Calculated from
%       calibration data & HumanBat_calibrate_cortex_ciholas
%
%   Outputs
%   -------
%   alignedCiholasCortex
%       Aligned ciholas and cortex data

cortex_fs = cortex.AnalogFrameRate(1);
ciholas_fs = ciholas.CDPmtdata.Fs;

posK_sample_global_ts_usec = ciholas.tag_data_filt{1}(:,8);
posM_sample_global_ts_usec = ciholas.tag_data_filt{2}(:,8);

% TODO: Check that ciholas and cortex are not missing samples / fill them
% in!

% Transform ciholas to cortex coordinates system
kq_pos = ciholas.tag_data_filt{1}(:,3:5);
ms_pos = ciholas.tag_data_filt{2}(:,3:5);
posK = ciholas2cortex.*kq_pos;
posM = ciholas2cortex.*ms_pos;
posK = resample(posK, cortex_fs, ciholas_fs); % resample to cortex 120hz
posM = resample(posM, cortex_fs, ciholas_fs); % resample to cortex 120hz
posK_resampled_global_ts_usec = resample(posK_sample_global_ts_usec,cortex_fs,ciholas_fs)*cortex_fs/ciholas_fs; % Global timestamps of resampled ciholas samples
posM_resampled_global_ts_usec = resample(posM_sample_global_ts_usec,cortex_fs,ciholas_fs)*cortex_fs/ciholas_fs; % Global timestamps of resampled ciholas samples

% At this point, posK, posM, and associated timestamps are in the
% coordinate system of Cortex with the same sampling frequency

% Slice data from first TTL
[minVal, ciholas_K_1st_ttl_indx] = min(abs(posK_resampled_global_ts_usec));
[minVal, ciholas_M_1st_ttl_indx] = min(abs(posM_resampled_global_ts_usec));
[minVal, cortex_1st_ttl_indx] = min(abs(cortex.global_sample_ts_usec));
cortex_1st_ttl_indx = cortex_1st_ttl_indx + 1.05*cortex_fs; % 1.05 sec offset between cortex and ciholas

alignedCortex.Markers = cortex.Markers(cortex_1st_ttl_indx:end,:,:);
alignedCortex.VideoFrameRate = cortex.VideoFrameRate;
alignedCortex.AnalogSignals = cortex.AnalogSignals(cortex_1st_ttl_indx:end,:);
alignedCortex.AnalogFrameRate = cortex.AnalogFrameRate;
alignedCortex.avgMarkerPos = cortex.avgMarkerPos(cortex_1st_ttl_indx:end,:);
alignedCortex.global_sample_ts_usec = cortex.global_sample_ts_usec(cortex_1st_ttl_indx:end);

alignedCiholas.fs = cortex_fs;
alignedCiholas.kq.pos = posK(ciholas_K_1st_ttl_indx:end,:);
alignedCiholas.kq.global_sample_ts_usec = posK_resampled_global_ts_usec(ciholas_K_1st_ttl_indx:end);
alignedCiholas.ms.pos = posM(ciholas_M_1st_ttl_indx:end,:);
alignedCiholas.ms.global_sample_ts_usec = posM_resampled_global_ts_usec(ciholas_M_1st_ttl_indx:end);

alignedCiholasCortex.cortex = alignedCortex;
alignedCiholasCortex.ciholas = alignedCiholas;
end