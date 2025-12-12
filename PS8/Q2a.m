%% Part (a): Replace transmitter in BigIdeal.m with m6params.m
% and analyze the spectrum of received signal r

clear all; close all; clc;

%% ===== MESSAGE GENERATION =====
% Create a test message
m = '0123456789 I wish I were an Oscar Meyer wiener 56789';

%% ===== FRAME PARAMETERS =====
% Define frame structure parameters
frameParams.userDataLength = 1;           % Number of symbols per frame
frameParams.preamble = '';                % No preamble for ideal system
frameParams.chanCodingFlag = 0;           % 0 = no channel coding
frameParams.bitEncodingFlag = 0;          % 0 = 8-bit ASCII (letters2pam)

%% ===== RF PARAMETERS =====
% Define RF and sampling parameters (NOTE: use underscores in field names!)
rfParams.f_s = 100;                       % Sampling frequency (Hz)
rfParams.T_t = 1;                         % Symbol period (seconds)
rfParams.T_t_err = 0;                     % Symbol period error (0 for ideal)
rfParams.f_if = 20;                       % Intermediate frequency (Hz)
rfParams.f_if_err = 0;                    % IF frequency error (0 for ideal)
rfParams.phaseNoiseVariance = 0;          % Phase noise variance (0 for ideal)
rfParams.SRRCLength = 4;                  % SRRC pulse length in symbols
rfParams.SRRCrolloff = 0.3;               % SRRC rolloff factor (beta)

%% ===== CHANNEL PARAMETERS =====
% Define channel characteristics
chanParams.c1 = [1 0 0];                  % Initial channel impulse response
chanParams.c2 = [1 0 0];                  % Final channel impulse response (time-invariant if c1=c2)
chanParams.randomWalkVariance = 0;        % Channel random walk variance (0 for ideal)
chanParams.SNR = Inf;                     % Signal-to-noise ratio (Inf = no noise)

% Adjacent User 1 parameters
chanParams.adjacentUser1Power = -Inf;    % Power in dB (-Inf = disabled)
chanParams.adjacentUser1f_if = 0;        % IF frequency of adjacent user 1
chanParams.adjacentUser1Chan = [1 0 0];  % Channel for adjacent user 1

% Adjacent User 2 parameters
chanParams.adjacentUser2Power = -Inf;    % Power in dB (-Inf = disabled)
chanParams.adjacentUser2f_if = 0;        % IF frequency of adjacent user 2
chanParams.adjacentUser2Chan = [1 0 0];  % Channel for adjacent user 2

% Narrowband interference parameters
chanParams.NBIfreq = 0;                   % NBI frequency
chanParams.NBIPower = -Inf;               % NBI power in dB (-Inf = disabled)

%% ===== TRANSMITTER: GENERATE RECEIVED SIGNAL =====
fprintf('Generating transmitted signal using B3IG...\n');
[r, s] = BigTransmitter(m, frameParams, rfParams, chanParams);

fprintf('Message length: %d characters\n', length(m));
fprintf('Number of transmitted symbols: %d\n', length(s));
fprintf('Received signal length: %d samples\n', length(r));
fprintf('Effective oversampling factor M: %d\n', round(length(r)/length(s)));

%% ===== ANALYZE SPECTRUM OF RECEIVED SIGNAL r =====

% Time vector for plotting
Ts = 1/rfParams.f_s;                      % Sampling period
t = (0:length(r)-1) * Ts;                 % Time vector

% Compute FFT and frequency vector
N = length(r);
R = fft(r);
f = rfParams.f_s * (0:N-1) / N;           % Frequency vector (0 to fs)
f_centered = rfParams.f_s * (-(N-1)/2:(N-1)/2) / N; % Centered frequency vector

% Shift FFT for centered spectrum
R_centered = fftshift(R);

%% ===== PLOTTING =====

figure('Position', [100 100 1200 800]);

% Plot 1: Time domain signal
subplot(3,2,1);
plot(t, r, 'b', 'LineWidth', 1);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Received Signal r(t) - Time Domain');
grid on;
xlim([0 min(5, max(t))]);  % Show first 5 seconds

% Plot 2: Magnitude spectrum (one-sided)
subplot(3,2,2);
plot(f(1:N/2), abs(R(1:N/2)), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of r(t) - One-Sided');
grid on;
hold on;
% Mark the IF frequency
% xline(rfParams.f_if, '--g', 'LineWidth', 2, 'Label', 'fIF');
legend('|R(f)|', 'IF Frequency');

% Plot 3: Magnitude spectrum (centered, two-sided)
subplot(3,2,3);
plot(f_centered, abs(R_centered), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of r(t) - Two-Sided (Centered)');
grid on;
hold on;
% xline(rfParams.f_if, '--g', 'LineWidth', 2);
% xline(-rfParams.f_if, '--g', 'LineWidth', 2);
legend('|R(f)|', '±fIF');

% Plot 4: Power Spectral Density
subplot(3,2,4);
plot(f(1:N/2), 20*log10(abs(R(1:N/2))+eps), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of r(t)');
grid on;
hold on;
% xline(rfParams.f_if, '--g', 'LineWidth', 2, 'Label', 'fIF');

% Calculate and display bandwidth
Rmag = abs(R(1:N/2));
Rmax = max(Rmag);
threshold = 0.01 * Rmax;  % 1% of maximum
bw_indices = find(Rmag > threshold);
if ~isempty(bw_indices)
    bandwidth = f(bw_indices(end)) - f(bw_indices(1));
    fprintf('\nEstimated signal bandwidth: %.2f Hz\n', bandwidth);
end

% Plot 5: Spectrogram (time-frequency analysis)
subplot(3,2,5);
window = hamming(256);
noverlap = 128;
nfft = 512;
spectrogram(r, window, noverlap, nfft, rfParams.f_s, 'yaxis');
title('Spectrogram of r(t)');
colorbar;

% Plot 6: Phase spectrum
subplot(3,2,6);
phase = angle(R_centered);
plot(f_centered, phase, 'm', 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Spectrum of r(t)');
grid on;

%% ===== DETAILED SPECTRUM ANALYSIS =====

fprintf('\n===== SPECTRUM ANALYSIS OF RECEIVED SIGNAL r =====\n');
fprintf('Sampling frequency f_s: %.2f Hz\n', rfParams.f_s);
fprintf('Nyquist frequency: %.2f Hz\n', rfParams.f_s/2);
fprintf('Intermediate frequency f_if: %.2f Hz\n', rfParams.f_if);
fprintf('Symbol rate: %.2f symbols/sec\n', 1/rfParams.T_t);
fprintf('SRRC rolloff factor: %.2f\n', rfParams.SRRCrolloff);

% Theoretical bandwidth calculation
symbol_rate = 1/rfParams.T_t;  % 1 Hz for T_t=1
theoretical_bw = symbol_rate * (1 + rfParams.SRRCrolloff);
fprintf('\nTheoretical signal bandwidth (baseband): %.2f Hz\n', theoretical_bw);
fprintf('Theoretical passband bandwidth: %.2f Hz\n', 2*theoretical_bw);

% Find peak in spectrum
[peak_val, peak_idx] = max(abs(R(1:N/2)));
peak_freq = f(peak_idx);
fprintf('\nPeak in spectrum at: %.2f Hz (expected at f_if = %.2f Hz)\n', peak_freq, rfParams.f_if);

% Signal power calculation
signal_power = sum(r.^2) / length(r);
fprintf('Average signal power: %.6f\n', signal_power);

%% ===== DISCUSSION =====
fprintf('\n===== DISCUSSION OF SPECTRUM r =====\n');
fprintf('1. The spectrum shows the 4-PAM signal modulated to IF = %.2f Hz\n', rfParams.f_if);
fprintf('2. Two-sided spectrum shows symmetric sidebands at ±f_if\n');
fprintf('3. SRRC pulse shaping creates a smooth spectral shape with rolloff = %.2f\n', rfParams.SRRCrolloff);
fprintf('4. Signal bandwidth is determined by symbol rate (%.2f Hz) and rolloff factor\n', symbol_rate);
fprintf('5. For ideal channel (SNR=Inf, no impairments), spectrum is clean\n');
fprintf('6. Signal is centered at f_if with no DC component (suppressed carrier)\n');
fprintf('7. Sampling frequency (%.2f Hz) is sufficient to avoid aliasing\n', rfParams.f_s);

%% ===== SAVE RESULTS =====
% Save the parameters and signal for part (b)
save('part_a_results.mat', 'r', 's', 'm', 'frameParams', 'rfParams', 'chanParams');
fprintf('\nResults saved to part_a_results.mat for use in part (b)\n');
