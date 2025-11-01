% -- Setup for Simulation--
fs = 1000;  %sampling frequency
Ts=1/fs;
t = 0:Ts:4; % Time vector from 0 to 4 seconds


%% --Create the Signal
x = sinc(100*(t-2));

%---plot the spectrum
plotspec(x,Ts);
title('Numerical Spectrumof x(t)=sinc (100(t-2))');