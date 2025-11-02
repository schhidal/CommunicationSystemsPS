% --- MATLAB/Octave Script for Part (a): Visualizing Upconversion ---

clear; close all;

% --- 1. Define Signal Parameters ---
f1 = 2000;         % 2 kHz, bandwidth of our baseband signal
f2 = 1000000;      % 1 MHz (1,000,000 Hz), our carrier frequency

% --- 2. Define Simulation Parameters ---
% Highest frequency in x2(t) will be f2 + f1 = 1,002,000 Hz.
% Nyquist rate is 2 * (f2 + f1) = 2,004,000 Hz.
% We must sample faster than this. Let's choose 2.5 MHz.
fs = 2.5e6;        % Sampling frequency (2,500,000 Hz)
Ts = 1/fs;         % Sampling interval
N = 125000;        % Number of samples to simulate
t = (0:N-1) * Ts;  % Time vector (0 to 0.05 seconds)

disp('Parameters set. Sampling at 2.5 MHz.');

% --- 3. Create the Baseband Signal x1(t) ---
% We'll make a stand-in signal by filtering white noise.
disp('Creating baseband signal x1(t)...');
% Generate white noise
white_noise = randn(1, N);
% Design a 6th-order Butterworth LPF with cutoff at f1 = 2 kHz
[b, a] = butter(6, f1 / (fs/2));
% Filter the noise to create our bandlimited signal
x1 = filter(b, a, white_noise);

% Plot the spectrum of x1(t)
figure;
plotspec(x1, Ts);
title('x_1(t) - Baseband Signal');
% Zoom in the frequency plot to see the baseband
subplot(2,1,2);
xlim([-f1*3, f1*3]); % Look at the +/- 6 kHz range

% --- 4. Create the RF Passband Signal x2(t) ---
disp('Creating RF passband signal x2(t)...');
c_rf = cos(2*pi*f2*t);    % 1 MHz Carrier
x2 = x1 .* c_rf;          % Modulate!

% Plot the spectrum of x2(t)
figure;
plotspec(x2, Ts);
% --- MATLAB/Octave Script for Part (a): Visualizing Upconversion ---

clear; close all;

% --- 1. Define Signal Parameters ---
f1 = 2000;         % 2 kHz, bandwidth of our baseband signal
f2 = 1000000;      % 1 MHz (1,000,000 Hz), our carrier frequency

% --- 2. Define Simulation Parameters ---
% Highest frequency in x2(t) will be f2 + f1 = 1,002,000 Hz.
% Nyquist rate is 2 * (f2 + f1) = 2,004,000 Hz.
% We must sample faster than this. Let's choose 2.5 MHz.
fs = 2.5e6;        % Sampling frequency (2,500,000 Hz)
Ts = 1/fs;         % Sampling interval
N = 125000;        % Number of samples to simulate
t = (0:N-1) * Ts;  % Time vector (0 to 0.05 seconds)

disp('Parameters set. Sampling at 2.5 MHz.');

% --- 3. Create the Baseband Signal x1(t) ---
% We'll make a stand-in signal by filtering white noise.
disp('Creating baseband signal x1(t)...');
% Generate white noise
white_noise = randn(1, N);
% Design a 6th-order Butterworth LPF with cutoff at f1 = 2 kHz
[b, a] = butter(6, f1 / (fs/2));
% Filter the noise to create our bandlimited signal
x1 = filter(b, a, white_noise);

% Plot the spectrum of x1(t)
figure;
plotspec(x1, Ts);
title('x_1(t) - Baseband Signal');
% Zoom in the frequency plot to see the baseband
subplot(2,1,2);
xlim([-f1*3, f1*3]); % Look at the +/- 6 kHz range

% --- 4. Create the RF Passband Signal x2(t) ---
disp('Creating RF passband signal x2(t)...');
c_rf = cos(2*pi*f2*t);    % 1 MHz Carrier
x2 = x1 .* c_rf;          % Modulate!

% Plot the spectrum of x2(t)
figure;
plotspec(x2, Ts);
title('x_2(t) - RF Passband Signal');
% Zoom in the frequency plot to the passband
subplot(2,1,2);
xlim([f2-f1*3, f2+f1*3]); % Look at the positive side
% If you want to see both sides, you can use:
% xlim([-f2-f1*3, f2+f1*3]);

disp('Done.');title('x_2(t) - RF Passband Signal');
% Zoom in the frequency plot to the passband
subplot(2,1,2);
xlim([f2-f1*3, f2+f1*3]); % Look at the positive side
% If you want to see both sides, you can use:
% xlim([-f2-f1*3, f2+f1*3]);

disp('Done.');