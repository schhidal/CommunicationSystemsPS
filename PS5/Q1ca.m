% --- MATLAB/Octave Script for DFT: Part (c) for signal (a) ---

clear; clc; close all;

% --- 1. Set up parameters (from specsin2.m) ---
f=100; Ts=1/1000; time=10.0;     % freq, sampling interval, time
N=2^10;                          % size of analysis window

% --- 2. Create the signal x[n] from part (a) ---
% We need a vector t from n=1 to n=N
t_vec = Ts:Ts:(N*Ts); % This matches the t(1:N) used in specsin2.m
                      % It creates t = [0.001, 0.002, ..., 1.024]
                      
% This is our new signal, the complex exponential:
% w = exp(j * 2*pi*f*t)
w = exp(1j * 2*pi*f * t_vec);

% --- 3. Compute the DFT (same as specsin2.m) ---
ssf=(-N/2:N/2-1)/(Ts*N);         % frequency vector
fw=fft(w);                       % do DFT/FFT
fws=fftshift(fw);                % shift it for plotting

% --- 4. Plot the results ---
figure(1);
% Plotting a complex signal in 2D isn't as clear as the sin wave,
% so we will just plot the real part to show it's a wave.
plot(t_vec, real(w));
xlabel('seconds'); ylabel('amplitude (Real Part)');
title('(a) Time Domain: x[n] = e^{j 2\pi f t}');

figure(2);
plot(ssf, abs(fws));             % plot magnitude spectrum
xlabel('frequency (Hz)'); ylabel('magnitude |X[k]|');
title('(a) Frequency Domain: DFT of Complex Sine Wave');