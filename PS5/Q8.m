% --- MATLAB/Octave Script for DC Removal Filter ---

clear; clc; close all;

% --- 1. Set up parameters ---
f_dc = 0;   % The DC component (0 Hz) - STOP
f1 = 100; % The cosine wave (100 Hz) - PASS

% To sample this, our sampling frequency (fs) must be
% more than 2 * (highest frequency) = 2 * 100 = 200 Hz.
% Let's use 1000 Hz.
fs = 1000;
Ts = 1/fs;
t = 0:Ts:(2-Ts); % Create 2 seconds of time

disp('--- Step 1 & 2: Setup and Signal Creation ---');
disp(['Sampling at ', num2str(fs), ' Hz.']);

% --- 2. Create the input signal x(t) ---
x = 1 + cos(2*pi*f1*t);

disp('Signal x(t) = 1 + cos(2*pi*100*t) created.');

% --- 3. Plot the spectrum of the input signal x(t) ---
disp('Plotting "Before" spectrum...');
figure; % Create a new figure window
plotspec(x, Ts); % Call your function to plot the time and frequency

% Add a title to the whole figure
title('"Before" Signal: x(t) = 1 + cos(2\pi*100*t)');

% --- (Optional) Zoom in on the frequency plot ---
subplot(2,1,2);
xlim([-150, 150]); % Only show frequencies from -150 to 150 Hz
% --- 4. Design the DC Removal (High-Pass) Filter ---
disp('Designing the DC Removal (HPF)...');

% Frequencies are normalized by Nyquist (fs/2 = 500 Hz)
f_nyq = fs/2;
f_edges = [0, 10, 50, f_nyq] / f_nyq;  % F = [0, 0.02, 0.1, 1.0]
a_edges = [0, 0,  1,   1];          % A = [0, 0, 1, 1]

% N = filter order (length - 1). 
% A larger N gives a sharper, "better" filter.
N_order = 50; 

% Design the filter 'b' (the impulse response)
b = firpm(N_order, f_edges, a_edges);

disp('Filter "b" has been designed.');

% --- (Optional) Plot the filter's frequency response ---
% 'freqz' plots the shape of the filter we just designed
figure; % Create a new figure (Figure 2)
freqz(b, 1, 512, fs); % Plot with 512 points, scaled to fs
title('DC Removal (HPF) Frequency Response');

% --- 5. Filter the signal x(t) with our HPF 'b' ---
disp('Filtering the signal...');

% 'filter' is the command to run the signal 'x' through our FIR filter 'b'.
% (The '1' is for IIR filters, so we just leave it as 1).
% This is the same 'filter' command from 'waystofilt.m'.
y = filter(b, 1, x);

disp('Signal has been filtered. Output is y.');

% --- 6. Plot the spectrum of the output signal y(t) ---
disp('Plotting "After" spectrum...');

figure; % Create a THIRD figure window
plotspec(y, Ts); % Plot the time and frequency of the output 'y'

% Add a title to the whole figure
title('"After" Signal: y(t) = HPF{x(t)}');

% --- (Optional) Zoom in on the frequency plot ---
subplot(2,1,2);
xlim([-150, 150]); % Use the same zoom as the "before" plot