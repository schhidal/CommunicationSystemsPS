% --- MATLAB/Octave Script for Convolution: Part (c) ---
% --- Plot and compare the spectra from (a) and (b) ---

clear; clc; close all;

% --- 1. Set up parameters ---
% Let's set T = 1 second for simplicity
T = 1;

% To get a good plot, we need a high sampling rate (fs)
% and a long time vector (num_seconds)
fs = 100;           % 100 samples per second
Ts = 1/fs;
num_seconds = 10;   % Simulate for 10 seconds
N = fs * num_seconds; % Total number of samples (1000)
t = (0:N-1) * Ts;     % Time vector from 0 to 9.99 s

% Create the frequency vector for plotting (from -fs/2 to fs/2)
f = (-N/2 : N/2-1) * (fs/N);

disp('--- Part (c): Plotting Spectra ---');
disp(['T = 1 sec, fs = 100 Hz, N = 1000 samples']);

% --- 2. Create Signal (a): The Single Pulse ---
disp('Creating signal (a): x_p(t)...');
x_p = zeros(1, N);
T_samples = T / Ts; % T = 1 sec, so T_samples = 100
% Put the pulse right in the middle (e.g., from t=4.5 to 5.5s)
start_index = (num_seconds/2 - T/2) * fs;
end_index = (num_seconds/2 + T/2) * fs - 1;
x_p(start_index : end_index) = 1;

% --- 3. Create Signal (b): The Square Wave ---
disp('Creating signal (b): x_w(t)...');
x_w = zeros(1, N);
T_p_samples = 2 * T / Ts; % Period is 2T = 200 samples
pulse_on = ones(1, T_samples);
pulse_off = zeros(1, T_samples);
one_period = [pulse_on, pulse_off];
% Repeat this pattern 5 times to fill the 1000-sample vector
x_w = repmat(one_period, 1, num_seconds / (2*T));

% --- 4. Compute the DFT/FFT of both signals ---
disp('Computing FFTs...');
% We scale by Ts*N to approximate the continuous FT magnitude
X_p_f = fftshift(fft(x_p)) * Ts; 
X_w_f = fftshift(fft(x_w)) * Ts;

% --- 5. Plot the magnitudes ---
disp('Plotting...');
figure;
subplot(2,1,1);
plot(f, abs(X_p_f));
title('(a) Spectrum of Single Pulse (X_p(f))');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([-4, 4]); % Zoom in to see the first few lobes
grid on;

subplot(2,1,2);
% We use 'stem' for the square wave because its spectrum is discrete
stem(f, abs(X_w_f));
title('(b) Spectrum of Square Wave (X_w(f))');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([-4, 4]); % Use the same zoom
grid on;