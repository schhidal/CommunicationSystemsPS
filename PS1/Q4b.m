clear; clc; close all;

% ---Parameters for the initial, oversampled signal---
f0 = 100;  % Frequency of the CT sinusoid in Hz
Ts_initial = 0.001;   % Initial sampling period in seconds in Hz
time_range = 0.1; % Time duration
t_initial = 0:Ts_initial:time_range-Ts_initial; % Initial time vector

% Generate the initial discrete-time signal x[n
x_initial = cos(2*pi*f0*t_initial);

% ---Plot the original (non-downsampled)spectrum for reference---
figure('Name', 'Initial Signal Spectrum');
plotspec(x_initial, Ts_initial);
title('Spectrum of Initial Sampled Signal (f_s = 1000 Hz)');
% Downsampling Factors
M_vals = [2, 3, 4, 5, 6];

% ---Loop, Downsample, and Plot ---
for i = 1:length(M_vals)
    M = M_vals(i);  % Current downsampling factor

    % Perform downsampling by keeping every M-th sample
    x_downsampled = x_initial(1:M:end);

    % The new sampling period is M times the original
    Ts_new = Ts_initial * M;
    fs_new = 1 / Ts_new;
figure;
    plotspec(x_downsampled, Ts_new);

    title_str = sprintf('Downsampled by M=%d (New f_s = %.1f Hz)', M, fs_new);
    title(title_str);
end

sgtitle('Spectra After Downsampling', 'FontSize', 14, 'FontWeight', 'bold');