% --- MATLAB/Octave Script for Convolution: Part (a) ---

clear; clc; close all;

% --- 1. Set up parameters ---
fs = 100;           % Sampling frequency (Hz)
Ts = 1/fs;          % Sampling interval
t = 0:Ts:(1-Ts);    % Time vector for 1 second (N=100 samples)
N = length(t);

% --- 2. Create the signals x[n] and y[n] ---
x = cos(2*pi*10*t); % Our 10 Hz cosine wave

y = zeros(1, N);    % A vector of zeros
y(1) = 1;           % A discrete-time impulse at n=0

% --- 3. Convolve in the Frequency Domain (using DFT) ---
% Length for linear convolution
N_conv = N + N - 1; % 100 + 100 - 1 = 199

% 'fft(signal, N_new)' automatically zero-pads the signal to length N_new
X = fft(x, N_conv); % X[k]
Y = fft(y, N_conv); % Y[k]

% Point-wise multiplication in frequency domain
Z = X .* Y;         % Z[k] = X[k] * Y[k]

% Inverse transform to get the result in time domain
z = ifft(Z);

% --- 4. Plot the results ---
figure;
% Time vector for the convolved signal (199 samples)
t_conv = (0:N_conv-1) * Ts;

plot(t_conv, z, 'b-', 'DisplayName', 'Result z(t)');
hold on;
plot(t, x, 'r--', 'DisplayName', 'Original x(t)');
title('(a) Convolution: cos(2\pi 10t) * \delta(t)');
xlabel('time (seconds)');
ylabel('amplitude');
legend;