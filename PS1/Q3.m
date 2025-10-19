% Low-Pass Filter
% Setuo Noise
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





