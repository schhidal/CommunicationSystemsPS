% --- MATLAB/Octave Script for DFT: Part (d) ---
% --- Shows the "leaky" vs. the "fixed" DFT ---

clear; clc; close all;

% --- 1. Set up parameters (from specsin2.m) ---
Ts=1/1000; time=10.0;
N=2^10;
t=Ts:Ts:(N*Ts); % Use only the N samples for the FFT
ssf=(-N/2:N/2-1)/(Ts*N);

% --- 2. "Leaky" Signal (from part c) ---
f_leaky = 100; % 100 Hz is NOT on a bin
w_leaky = cos(2*pi*f_leaky*t);
fws_leaky = fftshift(fft(w_leaky));

% --- 3. "Fixed" Signal (Our new magic frequency) ---
% Find the frequency for bin m=100
m = 100;
delta_f = 1 / (Ts * N); % This is 1/1.024
f_fixed = m * delta_f; % This is 100 / 1.024 = 97.65625 Hz

w_fixed = cos(2*pi*f_fixed*t);
fws_fixed = fftshift(fft(w_fixed));

% --- 4. Plot the comparison ---
figure;
subplot(2,1,1);
plot(ssf, abs(fws_leaky));
title(['(c) "Leaky" DFT (f = 100 Hz, 102.4 cycles)']);
ylabel('magnitude');
xlim([-150, 150]); % Zoom in

subplot(2,1,2);
plot(ssf, abs(fws_fixed));
title(['(d) "Fixed" DFT (f = 97.656 Hz, 100 cycles)']);
xlabel('frequency (Hz)'); ylabel('magnitude');
xlim([-150, 150]); % Zoom in