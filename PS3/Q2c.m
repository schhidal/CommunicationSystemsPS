% --2(b) Super-Nyquist Sampling--
fc = 1000; % carrier frequency
fs = (5/6) * fc; % oversampling fs = 10*fc = 10000 Hz
Ts = 1/fs; % Sampling interval in seconds
t = 0:Ts:4;  % Time vector from 0 to 4 seconds

disp('Parameters set. Sampling at 10,000 Hz.');

% --Create the RF signal y(t)--
disp('Generating RF signal y(t)...');
x_basebad = sinc(100*(t-2)); % Baseband signal
c_rf = cos(2*pi*fc*t);  % 1000 Hz Carrier
y_rf = x_basebad .* c_rf; % Modulated RF signal

% -- Plot the spectrum using your plotspec.m --
disp('Plotting spectrum....');
figure;
plotspec(y_rf, Ts);

% Manually adjust the title on the bottom plot for clarity
subplot(2,1,2);
title('(b) Super-Nyquist Spectrum (fs = 10 Khz)');
% Set x-axis limits to clearly see the components
xlim([-fs/2, fs/2]);



%% This Code is for one by one testing. Just required to change the Sampling frequency fs