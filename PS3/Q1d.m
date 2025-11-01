% -- Define signal parameters
fc = 1000;    % Carrier frequency in Hz
f_if = 200;    % Target INtermediate frequency  in Hz
f_baseband_max = 50;    % Max frequency of our Baseband sinc signal in Hz

% -- Define Simulation Parameters --
fs = 5000; % Samling frequency in Hz
Ts = 1/fs; % Sampling Time interval
t = 0:Ts:4; % Time vector from 0 to 4 seconds

disp('Parameters set. Sampling at 5000 Hz.');

% -- Create the original signal RF from y(t) part 1c
disp('Generating RF signal y(t)...');
x_baseband = sinc(100*(t-2));  % Baseband signal
c_rf = cos (2*pi*fc*t); % 1000 Hz Carrier
y_rf = x_baseband .* c_rf;  %  Modulated RF Signal

% -- Perform Low side injection Downconversion --
disp('Performing Low-side injection...');
f_lo_low = fc - f_if; % 1000-800 =200 Hz
c_lo = cos(2*pi*f_lo_low*t);  % 800 Hz local oscillator
w_lo = y_rf .* c_lo;  % The new downconverted signal

% --Plot the result using your plotspec.m --
figure;
plotspec(w_lo, Ts);

% --Manually adjust the title on the bottom for clarity --
subplot(2,1,2);
title('1(d) Spectrum after Low-side Injection (F\LO = 800 Hz)');
%-- Set x-axis limits to clearly all components
xlim([-fs/2, fs/2]);


