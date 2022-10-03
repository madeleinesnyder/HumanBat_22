% Running the script  hippo_extract_csc.m = Extract CSC data (EEG) and save into mat-file.
% See details in hippo_extract_csc.m .


clear all ; close all ; fclose all ; pack ;

%-----------------
% Michael Yartsev
%-----------------

% ##############################
% ##### Bat 6053 (Glubshi) #####
% ##############################


% ============= Parameters -- for a particular Day-Tetrode combination: ==========================
% ======= Day 11 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day11_neural_data_and_calibration\2009-05-20_09-35-18_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day11\VT_extracted_bat6053_Day11.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 39935921  1026276819 ]; %This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1211245885 2415753442 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2614238614  3614222220 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2280 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day11\CSC_extracted_bat6053_Day11_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

% ======= Day 12 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day12_neural_data_and_calibration\2009-05-21_07-20-41_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day12\VT_extracted_bat6053_Day12.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 1657843003  2598342472 ]; % This is currently taken from the event timestmaps during recording
%We take off 1 Min here from the behav session - see event strings 
times_microsec_behav_session_2 = [ 2863737472.00000,4333640357.00000 - 60*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 4668555296  5605945602 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2320 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day12\CSC_extracted_bat6053_Day12_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------


% ======= Day 13 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day13_neural_data_and_calibration\2009-05-22_12-04-41_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day13\VT_extracted_bat6053_Day13.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 251401394  1194763121 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1354412835.00000,3192335667.00000 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 3627657652  4638471810 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2380 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day13\CSC_extracted_bat6053_Day13_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------


% ======= Day 14 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day14_neural_data_and_calibration\2009-05-23_08-54-51_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day14\VT_extracted_bat6053_Day14.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 829677196  1748944339 ]; % This is currently taken from the event timestmaps during recording
%We take off 2 Min (120 sec) here from the behav session - see event strings 
times_microsec_behav_session_2 = [ 1961740391.00000,3763206807.00000 - 120*10^6  ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 3783149332 4765629688 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2400 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day14\CSC_extracted_bat6053_Day14_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------


% ======= Day 15 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day15_neural_data_and_calibration\2009-05-24_11-46-41_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day15\VT_extracted_bat6053_Day15.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 1928402763,2887494625 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 3159545127.00000,4703129625.00000 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 4866532142,6152203593 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [3257995175 3292337523]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2240 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day15\CSC_extracted_bat6053_Day15_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------


% ======= Day 16 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day16_neural_data_and_calibration\2009-05-25_10-10-02_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day16\VT_extracted_bat6053_Day16.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 3436903640,4359237340 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 4512635433,5923076706 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 6117789904,7079689492 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2240 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day16\CSC_extracted_bat6053_Day16_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------


% ======= Day 17 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day17_neural_data_and_calibration\2009-05-26_08-44-04_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day17\VT_extracted_bat6053_Day17.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 48288245,942749791 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1146883077,2451195842 ]; % This is currently taken from the event timestmaps during recording
%I am taking only the last 15 minutes of the post sleep session - see
%comments in the word file of today's recordings
times_microsec_post_sleep_session_3 = [ 4511766580-15*60*10^6  4511766580]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2320 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day17\CSC_extracted_bat6053_Day17_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------


% ======= Day 21 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day21_neural_data_and_calibration\2009-05-29_07-17-04_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day21\VT_extracted_bat6053_Day21.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 31150154,731794008 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 866239774,2140979239 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2704990708  3720087737]; % This is currently taken from the event timestmaps during recording
% All of the data we throw away below are from noises we had in EEG today
% during the POST SLEEP SESSION. See event file for details.
t_throw_away_data = [3147073531-30*10^6 3147073531; 3167564900-30*10^6 3167564900; ...
    3190949702-30*10^6 3190949702; 3251128937-20*10^6 3251128937; 327020239730-20*10^6 3270202397; ...
    3619309545-30*10^6 3619309545]; 
tetrode_depth_nominal_microns = 2380 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day21\CSC_extracted_bat6053_Day21_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

% ======= Day 22 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day22_neural_data_and_calibration\2009-05-30_06-32-21_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day22\VT_extracted_bat6053_Day22.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 56040127,891085351 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1012979820,2450038560-30*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [2634186346,3657584121 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2380 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day22\CSC_extracted_bat6053_Day22_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------
clear all ; close all ; fclose all ; pack ;

% ======= Day 23 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day23_neural_data_and_calibration\2009-05-31_10-36-12_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day23\VT_extracted_bat6053_Day23.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 137450251,1616110416 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1692310552,2521537937 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [2644359180,3973138736 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2400 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day23\CSC_extracted_bat6053_Day23_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;

% ======= Day 24 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day24_neural_data_and_calibration\2009-06-01_08-17-18_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day24\VT_extracted_bat6053_Day24.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 761230192,1715168898 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1842253577,3065914424-60*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 3408920877,4311643724 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2400 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day24\CSC_extracted_bat6053_Day24_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 25 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day25_neural_data_and_calibration\2009-06-02_06-09-37_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day25\VT_extracted_bat6053_Day25.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 41498568,910184149 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1051469453,2069769716 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2292547748,3310916765 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2420 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day25\CSC_extracted_bat6053_Day25_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------
clear all ; close all ; fclose all ; pack ;

% ======= Day 26 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day26_neural_data_and_calibration\2009-06-03_09-39-28_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day26\VT_extracted_bat6053_Day26.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 21076463,935431545 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1221382738,2332443729-60*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2521698676,3432079661 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2420 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day26\CSC_extracted_bat6053_Day26_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------
clear all ; close all ; fclose all ; pack ;


% ======= Day 27 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day27_neural_data_and_calibration\2009-06-04_10-52-04_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day27\VT_extracted_bat6053_Day27.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 64395226,997013288 - 60*10^6  ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1279890382,2252656098 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2414433623,3439429775 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2420 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day27\CSC_extracted_bat6053_Day27_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------
clear all ; close all ; fclose all ; pack ;

% % ======= Day 28 - For some reason the Low-Cut (the High-Pass Filter) was set to 10Hz today, so we can't use it for theta, hence I will not analyze this day ==================
% 
% filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day28_neural_data_and_calibration\2009-06-05_09-33-01_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
% filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day28\VT_extracted_bat6053_Day28.mat' ; % VT file: I will extract some data from it
% times_microsec_pre_sleep_session_1 = [ 381462186,1323053905  ]; % This is currently taken from the event timestmaps during recording
% times_microsec_behav_session_2 = [ 1564516786,2643653938 ]; % This is currently taken from the event timestmaps during recording
% times_microsec_post_sleep_session_3 = [ 2846611103,3828683747 ]; % This is currently taken from the event timestmaps during recording
% t_throw_away_data = []; 
% tetrode_depth_nominal_microns = 2420 ; % Copied from the file "comments_Glubshi_6053.doc"
% filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day28\CSC_extracted_bat6053_Day28_tt3.mat' ;
% 
% % ---------------
% hippo_extract_csc_Michael ; % Run the program itself
% % ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 29 ==================

filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat6053_behav_neural\Day29_neural_data_and_calibration\2009-06-07_08-38-49_neural_data_and_calibration\EEG3.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day29\VT_extracted_bat6053_Day29.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 26830545,1035414133  ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1197830341,2147548268 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2364747360,3302141572 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; 
tetrode_depth_nominal_microns = 2440 ; % Copied from the file "comments_Glubshi_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat6053_behav_neural\Day29\CSC_extracted_bat6053_Day29_tt3.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ##############################
% ##### Bat 7545 (No-Name) #####
% ##############################


% ======= Day 13 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day13_neural_data_and_calibration\2009-05-27_11-57-14_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day13\VT_extracted_bat7545_Day13.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 522407591  1159437567 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1823117960 3658738042 - 60*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 3807502518  4744270103 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [2572881631-60*10^6 2572881631;2842207641-60*10^6 2842207641;3139461229-60*10^6 3139461229]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2860 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day13\CSC_extracted_bat7545_Day13_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 14 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day14_neural_data_and_calibration\2009-05-28_08-03-39_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day14\VT_extracted_bat7545_Day14.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 1275000080, 2195550888 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 2324825944 3521015950 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 3682548479,4615870190 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = []; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2880 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day14\CSC_extracted_bat7545_Day14_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 15 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day15_neural_data_and_calibration\2009-05-29_08-38-32_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day15\VT_extracted_bat7545_Day15.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 65827492,997852439 ]; % This is currently taken from the event timestmaps during recording
%We remove the last 30 sec of behav which took until I wrote down the end
%of behav comment in the events.
times_microsec_behav_session_2 = [ 1095979353 2505964481 - 30*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2645113163,3553268429 ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [1587243701-60*10^6 1587243701; 1809206422-60*10^6 1809206422; ...
    1841828153-60*10^6 1841828153]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2900 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day15\CSC_extracted_bat7545_Day15_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;

% ======= Day 16 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day16_neural_data_and_calibration\2009-05-30_08-20-49_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day16\VT_extracted_bat7545_Day16.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 465491917,1370426022 ]; % This is currently taken from the event timestmaps during recording
%We remove the last 30 sec of behav which took until I wrote down the end
%of behav comment in the events.
times_microsec_behav_session_2 = [1516538261,2745028274-30*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2886588815,3806725412  ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [ ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2940 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day16\CSC_extracted_bat7545_Day16_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 17 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day17_neural_data_and_calibration\2009-05-31_14-09-02_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day17\VT_extracted_bat7545_Day17.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 207804263,934839619 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1041866345,1824130433 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2062831097,2607692452  ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [ ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2880 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day17\CSC_extracted_bat7545_Day17_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 20 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day20_neural_data_and_calibration\2009-06-03_08-16-58_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day20\VT_extracted_bat7545_Day20.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 37523705,955251424 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1065243097,2114345986 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2278904706,3249510691  ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [ ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2520 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day20\CSC_extracted_bat7545_Day20_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 21 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day21_neural_data_and_calibration\2009-06-04_08-02-39_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day21\VT_extracted_bat7545_Day21.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 45881521,943746576 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1037421464,2022497387 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2201150273,3160438001  ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [ ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2520 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day21\CSC_extracted_bat7545_Day21_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

clear all ; close all ; fclose all ; pack ;
% ======= Day 22 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day22_neural_data_and_calibration\2009-06-05_07-28-19_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day22\VT_extracted_bat7545_Day22.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 51959574,1047182810-60*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1394261568,2318378259-30*10^6 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2521417849,3509601238   ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [ ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2500 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day22\CSC_extracted_bat7545_Day22_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------
clear all ; close all ; fclose all ; pack ;

% ======= Day 26 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day26_neural_data_and_calibration\2009-06-09_08-00-53_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day26\VT_extracted_bat7545_Day26.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 172258205,1112743356 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1189094340 2663290137]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2836583371,3988226565   ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [ 1370271058, 1370271058-120*10^6; 1543436972,1667363620 ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2360 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day26\CSC_extracted_bat7545_Day26_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

% ======= Day 27 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day27_neural_data_and_calibration\2009-06-10_08-00-23_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day27\VT_extracted_bat7545_Day27.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 162100562,1075196482 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1368075287,2916281288 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 3054770809,4082511380  ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [  ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2360 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day27\CSC_extracted_bat7545_Day27_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

% ======= Day 28 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day28_neural_data_and_calibration\2009-06-11_08-00-10_neural_data_and_calibration\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day28\VT_extracted_bat7545_Day28.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 83875963,1048957006 ]; % This is currently taken from the event timestmaps during recording
times_microsec_behav_session_2 = [ 1222276132, 2633488380 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2771306997,3932310726  ]; % This is currently taken from the event timestmaps during recording
t_throw_away_data = [1717031447-180*10^6 1717031447 ]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2360 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day28\CSC_extracted_bat7545_Day28_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

% ======= Day 38 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day38_neural_data_and_live_video\2009-06-25_10-59-08_neural_data_and_live_video\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day38\VT_extracted_bat7545_Day38.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [  ]; % No Sleep session today
times_microsec_behav_session_2 = [ 231059735, 1108583517 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [   ]; % No Sleep session today
t_throw_away_data = [710623138 713100576]; % Throw away times when the bat DID NOT BEHAVE 
tetrode_depth_nominal_microns = 2240 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day38\CSC_extracted_bat7545_Day38_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael_NO_RIPPLES ; % Run the program itself
% ---------------

% ======= Day 39 ==================
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day39_neural_data_and_live_video\2009-06-26_09-23-32_neural_data_and_live_video\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day39\VT_extracted_bat7545_Day39.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 36195167,881312885 ]; % 
times_microsec_behav_session_2 = [ 1107586886 1887696770 ]; % This is currently taken from the event timestmaps during recording
times_microsec_post_sleep_session_3 = [ 2021553653,2853505957  ]; % 
t_throw_away_data = [  ]; % Throw away times when the bat DID NOT BEHAVE - This is based today on the live video we recorded
tetrode_depth_nominal_microns = 2240 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day39\CSC_extracted_bat7545_Day39_tt4.mat' ;

% ---------------
hippo_extract_csc_Michael ; % Run the program itself
% ---------------

% ======= Day 40 ==================
filename_CSC_AUDIO_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day40_neural_data_and_live_video\2009-07-08_12-38-03_neural_data_and_live_video\Audio1.Ncs'; % CSC-AUDIO (mic) file to extract
filename_CSC_EEG_in = 'D:\Michael\Data\Exp_Data\yr2009_bat7545_behav_neural\Day40_neural_data_and_live_video\2009-07-08_12-38-03_neural_data_and_live_video\EEG4.Ncs'; % CSC-EEG file to extract
filename_associated_VT_file = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day40\VT_extracted_bat7545_Day40.mat' ; % VT file: I will extract some data from it
times_microsec_pre_sleep_session_1 = [ 108851363,727102413-120*10^6 ]; % 
times_microsec_behav_session_2 = [ 998160682 2212694165]; % This is taken from the event live video file we took today
times_microsec_post_sleep_session_3 = [ 2386326755,2877135789  ]; % 
t_throw_away_data = [  ]; % Throw away times when the bat DID NOT BEHAVE - This is based today on the live video we recorded
tetrode_depth_nominal_microns = 2240 ; % Copied from the file "comments_noname_6053.doc"
filename_out = 'D:\Michael\Data\Expdata_Processed\yr2009_bat7545_behav_neural\Day40\CSC_extracted_bat7545_Day39_tt4.mat' ;

% ---------------
hippo_extract_EEG_and_audio_Michael ; % Run the program itself
% ---------------






% ======= THE END ===============

