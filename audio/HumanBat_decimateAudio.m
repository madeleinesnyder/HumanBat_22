function [outputArg1,outputArg2] = HumanBat_decimateAudio(concatenatedAudio)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clear ttlConCat;
    audioConCat_24khz = decimate(audioConCat,8);
    clear audioConCat;
    ttlConCat = load('ttlConCat_segment_0.mat');
    ttlConCat_24khz = decimate(ttlConCat,8);
    fs = 24000;
    mkdir([filePath 'processed']);
    save(fullfile(filePath,'processed/decimatedAudio.mat'), 'ttlConCat_24khz', 'audioConCat_24khz', 'fs', '-v7.3')
 
outputArg1 = inputArg1;
outputArg2 = inputArg2;
endfunction [outputArg1,outputArg2] = HumanBat_decimateAudio(concatenatedAudio)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clear ttlConCat;
    audioConCat_24khz = decimate(audioConCat,8);
    clear audioConCat;
    ttlConCat = load('ttlConCat_segment_0.mat');
    ttlConCat_24khz = decimate(ttlConCat,8);
    fs = 24000;
    mkdir([filePath 'processed']);
    save(fullfile(filePath,'processed/decimatedAudio.mat'), 'ttlConCat_24khz', 'audioConCat_24khz', 'fs', '-v7.3')
 
outputArg1 = inputArg1;
outputArg2 = inputArg2;