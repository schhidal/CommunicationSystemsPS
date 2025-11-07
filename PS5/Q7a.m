% --- MATLAB/Octave Script for Notch Filter: Part (a) ---

clear; clc; close all;

% --- 1. Set up parameters ---
f1 = 280; % Frequency 1 (280 Hz) - PASS
f0 = 300; % Frequency 2 (300 Hz) - STOP
f2 = 320; % Frequency 3 (320 Hz) - PASS

fs = 1000; % Sampling frequency (must be > 2*f2)
Ts = 1/fs;
t = 0:Ts:(2-Ts); % Create 2 seconds of time
disp('--- Part (a): Pole-Zero Placement Method ---');
disp(['Sampling at ', num2str(fs), ' Hz.']);

% --- 2. Create the input signal x(t) ---
x = cos(2*pi*f1*t) + cos(2*pi*f0*t) + cos(2*pi*f2*t);
disp('Signal x(t) created with 280, 300, and 320 Hz tones.');

% --- 3. Plot the "Before" spectrum of x(t) ---
figure; % Figure 1
plotspec(x, Ts);
title('"Before" Signal: x(t) with 3 tones');
subplot(2,1,2);
xlim([-fs/2, fs/2]); % Show positive frequencies

%--- 4. Design the notch fitlerf 
disp('Designing the 300 Hz Notch Filter...');
r=0.95; % Pole radius (closer to 1 = narror notch)
Omega_0 = 2*pi*f0/fs;    % Digital frequency to notch 

% 'b' coefficients (from zeros)
b_unscaled = [1, -2*cos(Omega_0), r^2];
% 'a' coefficients (from poles)
a = [1, -2*r*cos(Omega_0), r^2];

% Calculate gain 'g' to make gain at DC = 1
g = (1 - 2*r*cos(Omega_0) + r^2) / (2 - 2 * cos(Omega_0));
b=g* b_unscaled;

disp('Filter "b" and "a" have been designed. ');
% ---- %. Plot the filter's frequency response

figure;
freqz(b,a, 512, fs);
title('(a) Pole-Zero Notch FIlter Response (r=0.95)');

%---6. Filter the signal x(t) with our filter
y= filter(b,a,x);  %%% Use IIR filter syntax: filter (b,a,data)
disp('Signal has been filtered. Ouptput is y.');

figure;
plotspec(y,Ts);
title('Filtered signal');
subplot(2,1,2);
xlim([-fs/2,fs/2]);

