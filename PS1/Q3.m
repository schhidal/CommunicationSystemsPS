% Low-Pass Filter
clear; clc; close all;
% Setup Noise
Fs = 20000;     % Sampling frequencz in Hz
Ts = 1/Fs;      % Sampling Period
T = 1;          % Total time duration in seconds
N = T / Ts;     % Total number of samples
t = (0:N-1) * Ts; % Time vector

% Generate uniform random noise between -1 and 1
noise = 2 * rand(1, N)-1;

% Plot Input Spectrum
figure('Name', 'LPF Input');
plotspec(noise, Ts);
title('Input Noise Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Design LPF using firpm
N_taps = 101; % Number of taps (Order = N_taps -1)
filter_order = N_taps - 1;
nyquist_freq = Fs / 2;

fc_lp = 5000; % cutoff frequency in Hz
trans_width_lp = 50; % Transition band width in Hz

% Define frequency bands (normalized by nyquist frequency)
% [0, Passband_End, Stopband_Start, Nyquist]
f_lp_norm = [0, fc_lp-trans_width_lp/2, fc_lp+trans_width_lp/2, nyquist_freq] / nyquist_freq;

% Define desired amplitudes in those bands
% [Pass, Pass, Stop, Stop]
a_lp = [1, 1, 0, 0];

% Design the filter coeffiecients
b_lp = firpm(filter_order, f_lp_norm, a_lp);

% Apply filter
filtered_noise_lp = filter(b_lp, 1, noise);

% Plot Output specterum
figure('Name', 'LPF Output');
plotspec(filtered_noise_lp, Ts);
title('Output Noise Spectrum After LPF (fc = 5 KHz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');




%% ----------3b ------------------------
% --- Part (b): High-Pass Filter ---
% (Add this code section after the code for part (a) or run independently
% after generating the 'noise' signal and setting Fs, Ts etc.)

fc_hp = 1000; % Cutoff frequency in Hz
trans_width_hp = 500; % Transition band width in Hz
nyquist_freq = Fs / 2; % Defined in part (a) code
filter_order = 100;    % Defined in part (a) code

% Define frequency bands (normalized by Nyquist frequency)
% [0, Stopband_End, Passband_Start, Nyquist]
f_hp_norm = [0, fc_hp-trans_width_hp/2, fc_hp+trans_width_hp/2, nyquist_freq] / nyquist_freq;

% Define desired amplitudes in those bands
% [Stop, Stop, Pass, Pass]
a_hp = [0, 0, 1, 1];

% Design the filter coefficients
b_hp = firpm(filter_order, f_hp_norm, a_hp);

% Apply Filter
filtered_noise_hp = filter(b_hp, 1, noise);

% Plot Output Spectrum
figure('Name', 'HPF Output');
plotspec(filtered_noise_hp, Ts);
title('Output Noise Spectrum After HPF (fc = 1 kHz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');




%-------3c------------------
% --- Part (c): Band-Pass Filter ---
% (This code assumes 'noise', 'Fs', 'Ts' are still in the workspace)

fl_bp = 2000; % Lower cutoff frequency in Hz
fh_bp = 5000; % Higher cutoff frequency in Hz
trans_width_bp = 500; % Transition band width in Hz
nyquist_freq = Fs / 2;
filter_order = 100;

% Define frequency bands (normalized by Nyquist frequency)
% [0, stop1_end, pass_start, pass_end, stop2_start, nyquist]
f_bp_norm = [0, fl_bp-trans_width_bp/2, fl_bp+trans_width_bp/2, ...
             fh_bp-trans_width_bp/2, fh_bp+trans_width_bp/2, nyquist_freq] / nyquist_freq;

% Define desired amplitudes in those bands
% [Stop, Stop, Pass, Pass, Stop, Stop]
a_bp = [0, 0, 1, 1, 0, 0];

% Design the filter coefficients
b_bp = firpm(filter_order, f_bp_norm, a_bp);

% Apply Filter
filtered_noise_bp = filter(b_bp, 1, noise);

% Plot Output Spectrum
figure('Name', 'BPF Output');
plotspec(filtered_noise_bp, Ts);
title('Output Noise Spectrum After BPF (Passband: 2-5 kHz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');




