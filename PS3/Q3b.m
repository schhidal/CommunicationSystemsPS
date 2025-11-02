% --- MATLAB/Octave Script for Part (b): BPF Visualization ---

clear; close all;

% --- 1. Define Signal Frequencies ---
f1 = 2000;         % 2 kHz (Baseband bandwidth)
f2 = 1000000;      % 1 MHz (Carrier)
f3 = 990000;       % 990 kHz (Interferer 1)
f4 = 1005000;      % 1005 kHz (Interferer 2)

% --- 2. Define BPF Cutoff Frequencies (Our Choice) ---
f5 = 995000;       % 995 kHz (BPF lower cutoff)
f6 = 1003000;      % 1003 kHz (BPF upper cutoff)

% --- 3. Define Simulation Parameters ---
% Highest frequency is f4. Nyquist rate is 2*f4 = 2,010,000 Hz.
% We'll use fs = 2.5 MHz, same as before.
fs = 2.5e6;        % Sampling frequency (2,500,000 Hz)
Ts = 1/fs;         % Sampling interval
N = 125000;        % Number of samples
t = (0:N-1) * Ts;  % Time vector (0 to 0.05 seconds)

disp('Parameters set. Sampling at 2.5 MHz.');

% --- 4. Create the signal components ---
disp('Generating signals...');
% a) Create the desired signal x2(t)
white_noise = randn(1, N);
[b_lpf, a_lpf] = butter(6, f1 / (fs/2)); % LPF for baseband
x1 = filter(b_lpf, a_lpf, white_noise);
c_rf = cos(2*pi*f2*t);
x2 = x1 .* c_rf; % Our desired signal

% b) Create the interferers
interferer1 = 0.5 * cos(2*pi*f3*t); % Interferer at 990 kHz (amplitude 0.5)
interferer2 = 0.5 * cos(2*pi*f4*t); % Interferer at 1005 kHz (amplitude 0.5)

% c) Add them all together to create x3(t)
x3 = x2 + interferer1 + interferer2;

% --- 5. Plot the spectrum of the "dirty" signal x3(t) ---
disp('Plotting "dirty" signal x3(t)...');
figure; % Create Figure 1
plotspec(x3, Ts);
title('x_3(t) - Signal BEFORE Filtering');
% Zoom in the frequency plot to the region of interest
subplot(2,1,2);
xlim([f3-1000, f4+1000]); % Look from 989 kHz to 1006 kHz

% --- 6. Design and apply the BPF to get x4(t) ---
disp('Designing and applying BPF...');
% Design a 6th-order Butterworth BPF
[b_bpf, a_bpf] = butter(6, [f5, f6] / (fs/2), 'bandpass');
x4 = filter(b_bpf, a_bpf, x3);

% --- 7. Plot the spectrum of the "clean" signal x4(t) ---
disp('Plotting "clean" signal x4(t)...');
figure; % Create Figure 2
plotspec(x4, Ts);
title('x_4(t) - Signal AFTER BPF');
% Use the same zoom settings
subplot(2,1,2);
xlim([f3-1000, f4+1000]);

disp('Done.');