% --- MATLAB / Octave Script for Low-Pass Filter Design

clear; clc; close all;

% -- Set up the configuration --

f1 = 100; % Frequency of the first signal
f2 = 300; % Frequency of the second signal

% To sample this, our sampling frequency fs must be more than 2 * highest frequency
% 2 * 300 = 600, Let's use 1000 Hz
fs = 1000
Ts = 1/fs;
t = 0:Ts:(2-Ts); % Create 2 seconds of time

disp('--- Step 1 & 2: Setup and signal Creation');
disp(['Sampling at ', num2str(fs), 'Hz. ']);


% -- 2. Create the input signal x(t) --
x = cos(2*pi*f1*t) + cos(2*pi*f2*t);

% -- 3. Plot the spectrum of the input signal x(t) --
figure;
plotspec(x,Ts);
title('Signal before filter: x(t) = cos(2\pi*100*t) + cos(2\pi*300*t)');
xlim([-fs/2, fs/2]);

% --4. Design the Low-Pass filter --
disp('Designing the 200 Hz LPF ...');

% Frequeniceis are normalized by Nyquist rate (fs/2 = 500 Hz)
f_edge = [1,150,250,500] / (fs/2); % F = [0, 0.3, 0.5, 1.0]
a_edge =[1,1,0,0]; % A = [1,1,0,0]

% N = filter order(length - 1)
% A larger N gives a sharper, better filter
N_order = 50,

% Design the filter 'b' (the impulse response)
b = firpm(N_order, f_edge, a_edge);

disp('Filter "b" has been designed. ');

% -- Plot the fileter's frequency response --
% 'freqz' plots the shape of the filter we just designed
figure;
freqz(b,1,512,fs); % Plot with 512 pints, points scaled to fs
title('LPF Frequency response');

% --5. FIlter the Signal --
disp('Filtering the signal ...');

% 'filter' is the command to run the signal 'x' through our FIR filter 'b'.
%  (The '1' is for IIR filters, so we just leave it as 1)
% This is the same 'filter' command from 'waystofilr.m' [cite: user request].
y = filter(b,1,x);

disp('Signal has been filtered. Output is y.');

% --6. Plot the spectrum of the output signal y(t) --
disp('Plotting spectrum of "after" signal y(t) ..... ');

figure; % Create a THIRD figure window
plotspec(y,Ts); % Plot the time and  frequency of the output 'y'

% Add a title to the whole figure
title('"After Signal : y(t) = LPF{x(t)}"');

% --(Optional) Zoom in on the frequency plot --
subplot(2,1,2);
xlim([-400, 400]); % Use the same zoom as the "before" plot







