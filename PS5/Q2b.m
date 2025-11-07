% --- MATLAB/Octave Script for Convolution: Part (b) ---

clear; clc; close all;

% --- 1. Set up parameters ---
fs = 100;           % Sampling frequency (Hz)
Ts = 1/fs;          % Sampling interval
t = 0:Ts:(1-Ts);    % Time vector for 1 second (N=100 samples)
N = length(t);

% --- 2. Create the signals x[n] and y[n] ---
x = cos(2*pi*10*t); % Our 10 Hz cosine wave

% Create the rectangular pulse y(t) = u(t) - u(t-0.1)
pulse_width_sec = 0.1;
pulse_width_samples = pulse_width_sec / Ts; % 0.1 / 0.01 = 10 samples
y = zeros(1, N);
y(1:pulse_width_samples) = 1;  % [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, ..., 0]

% --- 3. Convolve in the Frequency Domain (using DFT) ---
% Length for linear convolution
N_conv = N + N - 1; % 199

X = fft(x, N_conv); % X[k]
Y = fft(y, N_conv); % Y[k]

% Point-wise multiplication in frequency domain
Z = X .* Y;         % Z[k] = X[k] * Y[k]

% Inverse transform to get the result in time domain
z = ifft(Z);
% NOTE: The result z will be complex, but the imaginary part
% should be tiny (e.g., 1e-15) and is just a numerical error.
% We take the real part.
z = real(z);

% --- 4. Plot the results ---
figure;
% Time vector for the convolved signal (199 samples)
t_conv = (0:N_conv-1) * Ts;

plot(t_conv, z, 'b-', 'DisplayName', 'Result z(t)');
hold on;
plot(t, x, 'r--', 'DisplayName', 'Original x(t)');
title('(b) Convolution: cos(2\pi 10t) * (u(t) - u(t-0.1))');
xlabel('time (seconds)');
ylabel('amplitude');
legend;
ylim([-2, 2]); % Set y-axis limits to compare