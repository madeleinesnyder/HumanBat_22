function [alignedCiholasCortexEphys,ephys_1st_ttl_indx] = HumanBat_alignCiholasCortexEphys(cortex, ciholas, ciholas2cortex, ephys_TTL, B_TT_unit)%,S_TT_unit)
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
ephys_fs = B_TT_unit(1).Fs;

% These time streams are the TTL?
posK_sample_global_ts_usec = ciholas.tag_data_filt{1}(:,8);
posM_sample_global_ts_usec = ciholas.tag_data_filt{2}(:,8);
posbat_sample_global_ts_usec = ciholas.tag_data_filt{7}(:,8);
accbat_sample_global_ts_usec = ciholas.tag_ac_data{7}(:,8);

% Shift ciholas - 1.05s
% Subtract 1.05 from the seconds colum (8) in the cortex tag tag data to
% aliugn the seconds in ciholas to the cortex



% TODO: Check that ciholas and cortex are not missing samples / fill them
% in!

% Transform ciholas to cortex coordinates system
kq_pos = ciholas.tag_data_filt{1}(:,3:5);
ms_pos = ciholas.tag_data_filt{2}(:,3:5);
bat_pos = ciholas.tag_data_filt{7}(:,3:5); 
%bat_acc1 = ciholas.tag_ac_data{7}(:,8);
bat_acc = ciholas.tag_ac_data{7}(:,3:5);
posK = ciholas2cortex.*kq_pos;              % align spatiallly
posM = ciholas2cortex.*ms_pos;              % align spatially
posbat = ciholas2cortex.*bat_pos;           % align spatially
accbat = ciholas2cortex.*bat_acc;           % align spatially
posK = resample(posK, cortex_fs, ciholas_fs); % resample ciholas data to cortex 120hz
posM = resample(posM, cortex_fs, ciholas_fs); % resample ciholas data to cortex 120hz
posbat = resample(posbat, cortex_fs, ciholas_fs); % resample ciholas data to cortex 120hz
accbat = resample(accbat, cortex_fs, ciholas_fs); % resample ciholas data to cortex 120hz
%accbat1 = resample(accbat1, cortex_fs, ciholas_fs); % resample ciholas data to cortex 120hz
posK_resampled_global_ts_usec = resample(posK_sample_global_ts_usec,cortex_fs,ciholas_fs)*cortex_fs/ciholas_fs; % Global timestamps of resampled ciholas samples
posM_resampled_global_ts_usec = resample(posM_sample_global_ts_usec,cortex_fs,ciholas_fs)*cortex_fs/ciholas_fs; % Global timestamps of resampled ciholas samples
posbat_resampled_global_ts_usec = resample(posbat_sample_global_ts_usec,cortex_fs,ciholas_fs)*cortex_fs/ciholas_fs; % Global timestamps of resampled ciholas samples
accbat_resampled_global_ts_usec = resample(accbat_sample_global_ts_usec,cortex_fs,ciholas_fs)*cortex_fs/ciholas_fs; % Global timestamps of resampled ciholas samples

% At this point, posK, posM, and associated timestamps are in the
% coordinate system of Cortex with the same sampling frequency

% Slice data from first TTL (first ciholas TTL is the absolute)
[minVal, ciholas_K_1st_ttl_indx] = min(abs(posK_resampled_global_ts_usec));
[minVal, ciholas_M_1st_ttl_indx] = min(abs(posM_resampled_global_ts_usec));
[minVal, ciholas_bat_1st_ttl_indx] = min(abs(posbat_resampled_global_ts_usec));
[minVal, ciholas_batacc_1st_ttl_indx] = min(abs(accbat_resampled_global_ts_usec));
[maxVal,ciholas_batacc_lst_ttl_indx] = max(abs(accbat_resampled_global_ts_usec));
[minVal, cortex_1st_ttl_indx] = min(abs(cortex.global_sample_ts_usec));
cortex_1st_ttl_indx = cortex_1st_ttl_indx + 1.05*cortex_fs; % 1.05 sec offset between cortex and ciholas
legnth of TTL vector / Fs of cortex = total amount of time between first and last ttl 
cortex_length = 

% Slice ephys data from the first TTL (first ciholas TTL is the absolute)
ttl_idxs = find(contains(ephys_TTL.event_types_and_details, 'Digital in. Digital input port status = 0x29ef0401; Digital input event status = 0x00000001; Digital in rising edge on pin number 1. Event name: Pin1Event;'));
ephys_1st_ttl_indx = ephys_TTL.event_timestamps_usec(ttl_idxs(1));
ephys_last_ttl_idx = ephys_TTL.event_timestamps_usec(ttl_idxs(end));

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
alignedCiholas.bat.pos = posbat(ciholas_bat_1st_ttl_indx:end,:);
alignedCiholas.bat.pos_global_sample_ts_usec = posbat_resampled_global_ts_usec(ciholas_bat_1st_ttl_indx:end);
alignedCiholas.bat.acc = accbat(ciholas_batacc_1st_ttl_indx:end,:);
alignedCiholas.bat.acc_global_sample_ts_usec = accbat_resampled_global_ts_usec(ciholas_batacc_1st_ttl_indx:end);
alignedCiholas.ttl.first = ciholas_batacc_1st_ttl_indx;
alignedCiholas.ttl.last = ciholas_batacc_lst_ttl_indx;


alignedEphysB.fs = ephys_fs;
for i=1:length(B_TT_unit)
    [val,idx] = min(abs(B_TT_unit(i).Timestamps-ephys_1st_ttl_indx));
    if (B_TT_unit(i).Timestamps(idx)-ephys_1st_ttl_indx) < 0
        idx = idx+1;
        val = B_TT_unit(i).Timestamps(idx);
    end
    shifted = B_TT_unit(i).Timestamps-ephys_1st_ttl_indx;
    shifted_chopped = shifted(idx:end);
    alignedEphysB.TT_unit_aligned(i).Timestamps = shifted_chopped./1e6;
    alignedEphysB.TT_unit_aligned(i).CellNumbers = B_TT_unit(i).CellNumbers(idx:end);    
    alignedEphysB.TT_unit_aligned(i).Features = B_TT_unit(i).Features(:,idx:end);
    alignedEphysB.TT_unit_aligned(i).NSpikes = length(alignedEphysB.TT_unit_aligned(i).Timestamps);
end

alignedCiholasCortexEphys.cortex = alignedCortex;
alignedCiholasCortexEphys.ciholas = alignedCiholas;
alignedCiholasCortexEphys.Bephys = alignedEphysB;
end