clear; clc; clear all;

% --- Define constants that do not change in the loop ---
f0 = 5; % Input signalfrequency in Hz (sin(2*pi*f0*t))
T =1; % Total time duration for simulation in seconds

% --- Define the carrier frequencies we want to test ---
fc_vals = [100, 1000, 10000]; % Vector of carrier frequencies in Hz

%% --- Loop through each carrier frequency ---
for i = 1:length(fc_vals)
    % Get the current carrier frequency for this iteration
    fc = fc_vals(i);

    % --- Step 1: Determine an adequate sampling frequency ---
    % To avoid aliasing, Fs must be > 2 * highest_frequency.
    % The highest frequency in our output will be fc + f0.
    % We'll use a safety factor of 4 for a clean plot.
    Fs = 4 * (fc + f0);
    Ts = 1/Fs;
    t = 0:Ts:T-Ts; % Time vector depends on Fs, so it's inside the loop

    % --- Step 2: Generate signals and mix them ---
    x = sin(2*pi*f0*t); % Input signal
    c = cos(2*pi*fc*t); % Carrier signal
    y = x .* c; % Modulated output signal

% --- Step 3: Plot the spectrum ---
figure;
plotspec(y, Ts);

% --- Step 4: Add a descriptive title and adjust axes ---
title_str = sprintf('Spectrum for f_c = %d Hz', fc);
title(title_str);
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Zoom in on the area of interest around the sidebands
xlim_range = fc + 200; % Set a viewing window around the carrier
xlim([-xlim_range, xlim_range]);
end

sgtitle('Analysis of a Mixer in the Frequency Domain', 'FontSize', 14, 'FontWeight', 'bold');