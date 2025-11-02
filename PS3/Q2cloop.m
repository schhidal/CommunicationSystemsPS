% --- MATLAB/Octave Script for Part (c): Sub-Nyquist Sampling Loop ---

clear; close all;
%-- 1.Define Signal and Simulation Parameters --
fc = 1000; % Carrier frequency in Hz
disp('Parameters set.');

% List of sampling frequency to Test
fs_fractions = [5/6, 3/4, 2/3, 1/2, 1/3, 1/4];
num_cases = length(fs_fractions);

% -- Lop through the each cases--
for i = 1:num_cases
    %---2a. Calculate new sampling parameters--
    fs = fs_fractions(i) * fc; % Sub-Nyquist sampliong frequency
    Ts = 1/fs;  % Sampling Interval in seconds

    t = 0:Ts:4;  % Time vector from 0 to 4

    fs_string = ['fs= ', num2str(fs_fractions(i), '%.2f'), '*fc = ', num2str(fs, '%.1f'), 'Hz'];
    disp(['Calculating for case ', num2str(i), ':',fs_string]);

    % -- 2b Create RF signal y(t) at the new sampling rate --
x_basebad = sinc(100*(t-2)); % Baseband signal
c_rf = cos(2*pi*fc*t);  % 1000 Hz Carrier
y_rf = x_basebad .* c_rf; % Modulated RF signal

% -- 2c Plot the spectrum using your plotspec.m --
figure;  % create a new figure for each case
plotspec(y_rf, Ts);

% Add a clear title to the whole figure
title(['(c) Sub-Nyquist Sampling with ', fs_string]);

% Adjust the bottom (frequency) plot's x-limits
subplot(2,1,2);
xlim([-fs/2, fs/2]);
end

disp('All 6 cases plotted.');
     