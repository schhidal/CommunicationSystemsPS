% --- Part (b): Simulate Ideal Suppressed Carrier AM ---

clear; clc; close all;

% --- Parameters ---
fm = 100;               % Message frequency in Hz
fc = 10000;             % Carrier frequency in Hz
Fs = 50000;             % Sampling frequency in Hz (Fs > 2*(fc+fm))
Ts = 1/Fs;              % Sampling interval
time = 0.1;             % Simulation duration in seconds
t = 0:Ts:time-Ts;       % Time vector (use 0 start for consistency)
N = length(t);          % Number of samples

% --- Transmitter ---
% Message signal
w = cos(2*pi*fm*t);
% Carrier signal
c = cos(2*pi*fc*t);
% Modulated signal (transmitted)
v = c .* w;

% --- Ideal Receiver ---
gam = 0; % Frequency offset (ideal = 0)
phi = 0; % Phase offset (ideal = 0)
% Receiver's local oscillator signal
c2 = cos(2*pi*(fc+gam)*t + phi);
% Signal after mixing at receiver
x = v .* c2;

% --- Low-Pass Filter Design ---
fl = 100;               % Filter order
nyquist_freq = Fs / 2;
% Pass frequencies up to slightly above fm (e.g., 200 Hz)
pass_end_norm = (fm + 100) / nyquist_freq;
% Stop frequencies starting well below 2*fc (e.g., 5 kHz)
stop_start_norm = 5000 / nyquist_freq;
fbe = [0, pass_end_norm, stop_start_norm, 1]; % Normalized band edges
damps = [1, 1, 0, 0];        % Desired amplitudes (Pass=1, Stop=0)
b_lpf = firpm(fl, fbe, damps);  % Design LPF coefficients

% --- Apply LPF to recover message ---
% The factor of 2 compensates for the 1/2 from mixing twice
m_recovered = 2 * filter(b_lpf, 1, x);

%% --- Plotting Spectra ---
figure('Name', 'Ideal AM Spectra');
ssf = (-N/2:N/2-1)*(Fs/N); % Frequency vector for plotting

subplot(4, 1, 1);
W_fft = fftshift(fft(w));
plot(ssf, abs(W_fft));
title('(1) Original Message Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([-fc*1.5 fc*1.5]); % Set consistent x-limits

subplot(4, 1, 2);
V_fft = fftshift(fft(v));
plot(ssf, abs(V_fft));
title('(2) Modulated (Transmitted) Signal Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([-fc*1.5 fc*1.5]);

subplot(4, 1, 3);
X_fft = fftshift(fft(x));
plot(ssf, abs(X_fft));
title('(3) Mixed Signal Spectrum (at Receiver, before LPF)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([-fc*1.5 fc*1.5]);

subplot(4, 1, 4);
M_rec_fft = fftshift(fft(m_recovered));
plot(ssf, abs(M_rec_fft));
title('(4) Recovered Message Spectrum (after LPF)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([-fc*1.5 fc*1.5]);

sgtitle('Ideal Suppressed Carrier AM (f_m=100 Hz, f_c=10 kHz)', 'FontSize', 14);

% Optional: Plot time-domain signals to visualize
% figure('Name', 'Ideal AM Time Domain');
% subplot(4,1,1); plot(t,w); title('Original Message'); xlim([0 5/fm]);
% subplot(4,1,2); plot(t,v); title('Modulated Signal'); xlim([0 5/fm]);
% subplot(4,1,3); plot(t,x); title('Mixed Signal (Receiver)'); xlim([0 5/fm]);
% subplot(4,1,4); plot(t,m_recovered); title('Recovered Message'); xlim([0 5/fm]);