function [pxx, f] = power_spectrum(data_dir, show_plots)
%check_ephys_noise Calculate power spectrum via welch's method of channel
%                  average csc
%   Parameters: data_dir | string | Directory containing extracted_data/
%               show_plots | boolean
data = load_extracted_data_logger(data_dir);

avg_csc = mean(data.csc, 2);
avg_csc = avg_csc(1:end-1);

fs = 31250; % Sampling frequency
f = round(fs/2);
f_res = 0.5; % Desired frequency resolution
T = 1.6/f_res; % Time window needed for desired freq res for hanning window;
window_size = fs*T; % Window size needed for desired time window
overlap_size = round(window_size / 2); % 50% window overlap
[pxx, f] = pwelch(avg_csc, window_size, overlap_size, fs/2, fs);

if(show_plots)
    figure;
    subplot(2,1,1);
    plot(avg_csc);
    xlabel('samples');
    ylabel('amplitude');
    title('Raw CSC trace');
    subplot(2,1,2);
    plot(f(1:100),pxx(1:100))
    set(gca, 'YScale', 'log')
    xlabel('Frequency (hz)');
    ylabel('Power (log)');
    title('Welchs Power Spectrum');
end
