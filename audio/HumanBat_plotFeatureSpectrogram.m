function HumanBat_plotFeatureSpectrogram(audio_clip) 

    % Plot all the SAP features using Ofer's toolbox from an audioclip

    load(audio_clip);
    short_seg = recbuf(:,1);
    
    fs = 192000;

    % Filter the raw signal without bandpass and see how the different filters look:
    % Smooth with a movmedian filter
    kk=20;
    MOV_mic_data = movmedian(short_seg,kk);
    %ABS_MOV_mic_data = downsample(tsmovavg(rms(abs(short_seg),2),'s',500,1),100);
    RMS_mic_data = tsmovavg(rms(abs(short_seg),2),'s',500,1); 
    figure('name',"Raw Data"); hold on; x1=subplot(4,1,1); hold on; title("Raw Data"); plot(short_seg); x2=subplot(4,1,2); hold on; title("MovMedian Filter"); plot(MOV_mic_data); x3=subplot(4,1,3); hold on; title("RMS"); plot(RMS_mic_data); linkaxes([x1 x2],'x');
    subplot(4,1,4); hold on; title("Periodogram")
    % Make periodograms
    periodogram(short_seg); hold off;    

    % From Markowitz using Ofer's software.
    % Extract peaks from m_AM
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]= sap_features(short_seg,fs);
    [pks,locs] = findpeaks(abs(m_AM),'MinPeakProminence',0.4);
    locs = locs*44;
    pkvect = NaN(1,length(short_seg));
    pkvect(locs) = 0.01;
    figure(); hold on; 
    subplot(4,1,1); hold on; plot(short_seg); plot(pkvect,'*r');
    subplot(4,1,2); hold on; plot(RMS_mic_data); plot(pkvect,'*r');
    subplot(4,1,4); hold on; 
    [IMAGE,F,T] = fb_pretty_sonogram(short_seg./abs(max(short_seg)),fs,'low',2.9,'zeropad',0);
    colormap(hot)
    imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
    ylabel('kHz')
    xlabel('time (s)');
    title('Spectrogram of manually selected audio subsegment');
    sound(short_seg*15,fs);

    %#-part spectrogram plot
    figure(); % make spectrogram..
    ax1 = subplot(3,1,1)

    [b,a]=ellip(5,.2,80,[100]/(fs/2),'high'); 
    %[b,a]=ellip(5,.2,80,[3515]/(fs/2),'high'); % filter above 3515
    [IMAGE,F,T] = fb_pretty_sonogram(filtfilt(b,a,short_seg./abs(max(short_seg))),fs,'low',2.9,'zeropad',0);
    colormap(hot)
    imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
    ylabel('kHz')
    xlabel('time (s)');
    title('Spectrogram of manually selected audio subsegment');

    % now plot RMS below
    [b,a]=butter(8,2*[40e3 80e3]./fs,'bandpass');
    convy=filter(b,a,short_seg);
    temp1 = zscore(zftftb_rms(convy',fs))*500;
    temp1(temp1<10) = 0; % get rid of filter artifacts..
    hold on;

    ax2 = subplot(3,1,2);
    plot((30:size(temp1,2))/fs,temp1(30:end),'b','LineWidth',3)
    %linkaxes([ax1,ax2],'x'); axis tight

    % plot the OG data
    ax3 = subplot(3,1,3); hold on; 
    plot((30:size(temp1,2))/fs,short_seg(30:end)); plot((30:size(temp1,2))/fs,pkvect(30:end),'*r');
    %plot(pk_vect,'*r'); 
    linkaxes([ax1,ax2,ax3],'x'); axis tight


end