% --- MATLAB/Octave Script for Convolution: Part (d) ---

clear; clc; close all;

% --- 1. Set up parameters ---
fs = 100;           % Sampling frequency (Hz)
Ts = 1/fs;          % Sampling interval
t = 0:Ts:(1-Ts);    % Time vector for 1 second (N=100 samples)
N = length(t);

% --- 2. Create the signals x[n] and y[n] ---
x = cos(2*pi*10*t); % Our 10 Hz cosine wave
y = sin(2*pi*10*t); % Our 10 Hz sine wave

% --- 3. Convolve in the Frequency Domain (using DFT) ---
% Length for linear convolution
N_conv = N + N - 1; % 100 + 100 - 1 = 199

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

% Create the theoretical answer: z(t) = 0.5 * t * sin(2*pi*10*t)
% We only care about the first 1 second (N samples)
t_theory = (0:N-1) * Ts;
% We create a 199-sample vector for comparison
z_theory_full = 0.5 .* t_conv .* sin(2*pi*10*t_conv);
% The theory is only valid while both signals overlap (first N samples)
% After N samples (1 sec), the convolution result changes.
% We will only plot the part where the theory is simple.

plot(t_conv, z, 'b-', 'DisplayName', 'Result z(t) from DFT');
hold on;
plot(t_theory, z_theory_full(1:N), 'r--', 'DisplayName', 'Theory: 0.5*t*sin(2\pi 10t)');
title('(d) Convolution: cos(2\pi 10t) * sin(2\pi 10t)');
xlabel('time (seconds)');
ylabel('amplitude');
legend;