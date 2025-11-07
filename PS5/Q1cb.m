% --- MATLAB/Octave Script for DFT: Part (c) for signal (b) ---
% --- This is specsin2.m, but with sin changed to cos ---

clear; clc; close all;

% --- 1. Set up parameters (from specsin2.m) ---
f=100; Ts=1/1000; time=10.0;     % freq, sampling interval, time
t=Ts:Ts:time;                    % define a time vector

% --- 2. Create the signal x[n] from part (b) ---
% THIS IS THE ONLY CHANGE
w=cos(2*pi*f*t);                 % define the sinusoid as a COSINE

% --- 3. Compute the DFT (same as specsin2.m) ---
N=2^10;                          % size of analysis window
ssf=(-N/2:N/2-1)/(Ts*N);         % frequency vector
fw=fft(w(1:N));                  % do DFT/FFT on the first N samples
fws=fftshift(fw);                % shift it for plotting

% --- 4. Plot the results (same as specsin2.m) ---
figure(1);
plot(t(1:N), w(1:N));             % plot the waveform
xlabel('seconds'); ylabel('amplitude');
title('(b) Time Domain: x[n] = cos(2\pi f t)');

figure(2);
plot(ssf, abs(fws));             % plot magnitude spectrum
xlabel('frequency (Hz)'); ylabel('magnitude |X[k]|');
title('(b) Frequency Domain: DFT of Cosine Wave');