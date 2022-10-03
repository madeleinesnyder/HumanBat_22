function HumanBat_plotSpectrogram(audio_clip,audio_loaded) 

    % Plot a pretty spectrogram of an audioclip .mat file from the audio
    % directory from B149f or B151
    
    if isempty(audio_clip)
        short_seg = audio_loaded;
    else
        load(audio_clip);
        short_seg = recbuf(:,1);
    end
    
    fs = 192000;
    if num2str(max(short_seg)) == '1'
        infinite_values = find(short_seg == max(short_seg));
        ii = NaN(size(short_seg,1),1);
        ii(infinite_values) = 0.5;
        %figure(); hold on; plot(short_seg); plot(ii,'*r');
        short_seg(infinite_values) = 0;
    end
    
    figure('name', 'Raw Data'); 
    [IMAGE,F,T] = fb_pretty_sonogram(short_seg./abs(max(short_seg)),fs,'low',2.9,'zeropad',0);
    colormap(hot)
    imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
    ylabel('kHz')
    xlabel('time (s)');
    title('Spectrogram of manually selected audio subsegment');

end