function [test_result] = test_ephys_noise(ephys_extracted_path, out_path)
%Test Ephys Noise
%   ephys_extracted_path : path to directory containing extracted_data/
%   folder from Julie's extract_logger_data function
%
%   out_path : path to test result output
[pxx, f] = power_spectrum(ephys_extracted_path, false);

fm = 600; % Max freq to plot

fig = figure();
tiledlayout(1,1);

ax1= nexttile;
plot(f(1:fm/2), 10*log10(pxx(1:fm/2)));
%title('');
xlabel('Frequency (Hz)');
ylabel('Power (dB/Hz)');
xticks([60 120 180 240 300 360 420 480 540 600]);
ylim([-20 50]);

sgtitle('Welchs Power Spectrum');
linkaxes([ax1], 'xy');
saveas(fig, fullfile(out_path));


test_result = true;

end

