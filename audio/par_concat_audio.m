fileNameStruct = dir('*_*_*.mat');
fname = fileNameStruct(1).name;
matches = regexp(fname, '_(.*)_*_.*.mat', 'tokens');

date = matches{1}{1};

fname_base = sprintf('_%s_audio_trial_',date);

chunkIdxConcat = [];
audioConcat = [];
parfor i=1:574
    fname = sprintf('%s%d.mat', fname_base, i);
    disp(fname)
    chunk = load(fname);
    chunk_num = str2num(chunk.audiofile(end));
    chunkIdxConcat = [chunkIdxConcat chunk_num];
    audioConcat = [audioConcat; chunk.recbuf(:,1)];
    ttlConcat = [ttlConcat; chunk.recbuf(:,end)];
end

audioConcat = load('concatenatedAudio.mat').audioConcat;
audioConcat = decimate(audioConcat,8);
save('mic1_38400hz.mat', "audioConcat");        

ttlConcat = load('concatenatedAudio.mat').ttlConcat;
ttlConcat = decimate(ttlConcat,8);

[R,LT,UT,LL,UL] = risetime(ttlConcat,24000);


ttlConcat = [];
ttlTimestamps = [];
ttlChunkNum = [];
for i=1:50
    fname = sprintf('%s%d.mat', fname_base, i);
    disp(fname)
    chunk = load(fname);
    chunk_num = str2num(chunk.audiofile(end));
    chunk_ttl = chunk.recbuf(:,end);
    [R,LT,UT,LL,UL] = risetime(chunk_ttl,chunk.fs);
    disp(size(LT))
    if(length(LT) <= 4)
        ttlTimestamps = [ttlTimestamps; LT];
        ttlChunkNum = [ttlChunkNum repmat([chunk_num], 1, length(LT))];
    end
    ttlConcat = [ttlConcat; chunk_ttl];
end
[R,LT,UT,LL,UL] = risetime(ttlConcat,chunk.fs);
