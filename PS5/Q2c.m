% --- MATLAB/Octave Script for Convolution: Part (c) ---

clear; clc; close all;

% --- 1. Set up parameters ---
fs = 100;           % Sampling frequency (Hz)
Ts = 1/fs;          % Sampling interval
t = 0:Ts:(1-Ts);    % Time vector for 1 second (N=100 samples)
N = length(t);

% --- 2. Create the signals x[n] and y[n] ---
x = cos(2*pi*10*t); % Our 10 Hz cosine wave

% Create the ramp pulse y(t) = 100*t for 0 <= t < 0.1
pulse_width_sec = 0.1;
pulse_width_samples = round(pulse_width_sec / Ts); % 10 samples
y = zeros(1, N);
% The time values for the ramp are [0, 0.01, ..., 0.09]
ramp_t = 0:Ts:(pulse_width_sec - Ts);
y(1:pulse_width_samples) = 100 * ramp_t;  % y = [0, 1, 2, ..., 9]

% --- 3. Convolve in the Frequency Domain (using DFT) ---
% Length for linear convolution
N_conv = N + N - 1; % 199

X = fft(x, N_conv); % X[k]
Y = fft(y, N_conv); % Y[k]

% Point-wise multiplication in frequency domain
Z = X .* Y;         % Z[k] = X[k] * Y[k]

% Inverse transform to get the result in time domain
z = ifft(Z);
z = real(z); % Discard tiny imaginary numerical errors

% --- 4. Plot the results ---
figure;
% Time vector for the convolved signal (199 samples)
t_conv = (0:N_conv-1) * Ts;

plot(t_conv, z, 'b-', 'DisplayName', 'Result z(t)');
hold on;
plot(t, x, 'r--', 'DisplayName', 'Original x(t)');
title('(c) Convolution: cos(2\pi 10t) * (ramp pulse)');
xlabel('time (seconds)');
ylabel('amplitude');
legend;