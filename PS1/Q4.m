clear; clc; clear all;
% -------Parameters-------
Ts = 0.001;   % Sampling period in seconds (fs = 1000 Hz)
time_range = 0.1;   % Represenatative time range in seconds
t = 0:Ts:time_range-Ts; % Time vector

f_vals = [400, 450, 500, 550, 600]; % Frequencies to test in Hz

%figure('Name', 'Part (a): Sampling and Aliasing');


% -----Loop and Plot-----
for i = 1:length(f_vals);
    f0 = f_vals(i);
    y = cos(2*pi*f0*t);  % Generate cosine wave
figure;
    %subplot(3,2,i); % Create a grid of Plots
    plotspec(y, Ts); % Use plotspec.m to plot time and frequency domain

    % Add and infromative title to each plot
    title_str = sprintf('Input f = %d Hz', f0);

    % Check for aliasing to add for the title
    nyquist_freq = 1 / (2*Ts);
    if f0 > nyquist_freq
        alias_freq = abs(f0 -round(f0/(1/Ts))*(1/Ts));
        title_str = [title_str, sprintf(' (Aliases to +/- %d Hz)', alias_freq)];
    elseif f0 == nyquist_freq
        title_str = [title_str,' (At Nyquist Frequency)'];
    end
    title(title_str);
end

sgtitle('Spectra of Sampled Cosine Waves (f_s = 1000 Hz)', 'Fontsize', 14, 'FontWeight', 'bold');