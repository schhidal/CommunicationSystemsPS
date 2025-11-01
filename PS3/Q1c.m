% -- Setup for simulation --
fc = 1000; % carrier frequency in Hz
fs = 4000; % New sampling frequency (must be >2*(1000+50)
Ts = 1/fs; % New sampling Interval
t = 0:Ts:4;  % Time vector from 0 to 4 Seconds

% --Create the BAseband Signal --
x = sinc(100*(t-2));

% -- Create the RF signal y(t) --
c = cos(2*pi*fc*t); % Carrier wave
y = x .* c; % Modulated RF Signal

% -- PLot the spectrum --
plotspec(y, Ts);
title('Numerical Spectrum of AM Signal y(t)');

