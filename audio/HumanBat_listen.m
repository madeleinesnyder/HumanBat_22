function HumanBat_listen(audioclip)

    % Plays raw audio clip ( .mat in audio folder from B149f or B151)
    % Saves manual input notes about contents of audio matfile to a .mat
    % for later curation (filename_e.mat indicates this clip has an
    % echolocation)

    load(audioclip);
    X = {};
    for i=1
        sound(recbuf(:,1),192000);
        prompt = 'Voice? Echolocations? Other? None? (V, E, O, N)';
        str = input(prompt,'s');
        X{i} = str;
    end
    save(strcat('__',audioclip(1:end-4),'_mic_notes_',X{1}),'X');
end