% --- MATLAB/Octave Script for Notch Filter: Part (c) ---

clear; clc; close all;

% --- 1. Set up parameters ---
f_pass1 = 280; % Frequency 1 (280 Hz) - PASS
f_stop = 300; % Frequency 2 (300 Hz) - STOP
f_pass2 = 320; % Frequency 3 (320 Hz) - PASS

fs = 1000; % Sampling frequency
Ts = 1/fs;
t = 0:Ts:(2-Ts); % Create 2 seconds of time
disp('--- Part (c): Chebyshev Type 1 Filter Method ---');
disp(['Sampling at ', num2str(fs), ' Hz.']);

% --- 2. Create the input signal x(t) ---
x = cos(2*pi*f_pass1*t) + cos(2*pi*f_stop*t) + cos(2*pi*f_pass2*t);
disp('Signal x(t) created with 280, 300, and 320 Hz tones.');

% --- 3. Plot the "Before" spectrum of x(t) ---
figure; % Figure 1
plotspec(x, Ts);
sgtitle('"Before" Signal: x(t) with 3 tones');
subplot(2,1,2);
xlim([0, fs/2]); % Show positive frequencies

% --- 4. Design the Notch Filter (using 'cheby1') ---
disp('Designing the 300 Hz Chebyshev-1 Notch Filter...');
f_nyq = fs/2;           % Nyquist frequency (500 Hz)
f_notch_low = 295;      % Low edge of our stop-band
f_notch_high = 305;     % High edge of our stop-band
N_order = 6;            % Filter order (same as before)
Rp = 1;                 % Passband ripple (1 dB)

% Normalized frequency band to stop
Wn = [f_notch_low, f_notch_high] / f_nyq; % [0.59, 0.61]

% Design the filter
[b, a] = cheby1(N_order, Rp, Wn, 'stop');

disp('Filter "b" and "a" have been designed.');

% --- 5. (Optional) Plot the filter's frequency response ---
figure; % Figure 2
freqz(b, a, 1024, fs);
title('(c) Chebyshev-1 Notch Filter Response (Order=6, Ripple=1dB)');

% --- 6. Filter the signal x(t) with our filter ---
disp('Filtering the signal...');
y = filter(b, a, x); % Use IIR filter syntax
disp('Signal has been filtered. Output is y.');

% --- 7. Plot the "After" spectrum of y(t) ---
disp('Plotting "After" spectrum...');
figure; % Figure 3
plotspec(y, Ts);
title('"After" Signal: y(t) = BPF{x(t)}');
subplot(2,1,2);
xlim([0, fs/2]); % Show positive frequencies