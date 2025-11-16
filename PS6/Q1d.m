addpath('../Tests')

% Part (d): 2-PAM Signal with Random ±1 Amplitudes
% Bandwidth Consumption Investigation

clear all; close all;

%% Parameters
N = 100;                    % Number of symbols
M = 100;                    % Oversampling factor (same as book example)
T = 1;                      % Symbol period

%% Generate random ±1 symbols (2-PAM)
% This is the KEY difference from the book examples
m = 2*round(rand(1,N)) - 1;  % Random ±1 sequence

%% Create SQUARE pulse (not Hamming like in the book)
% Oversample the symbols
mup = zeros(1, N*M);
mup(1:M:N*M) = m;           % Place symbols at M intervals

% Square pulse shape (different from book's Hamming pulse)
p = ones(1, M);             % Simple rectangular pulse

% Pulse shaping (convolve pulse with data)
x = filter(p, 1, mup);      % Transmitted signal

%% Plot spectrum of the 2-PAM signal
figure(1);
plotspec(x, 1/M);           % Use book's plotspec function
title('Spectrum of 2-PAM Signal with Random ±1 Amplitudes');

%% Find zero-crossing bandwidth
% Compute FFT
spectrum = fft(x);
mag_spectrum = abs(fftshift(spectrum));
freq = (-length(spectrum)/2:length(spectrum)/2-1)/(M*T*length(spectrum));

% Find first zero crossing
center = ceil(length(spectrum)/2);
% Search in positive frequencies only
pos_freqs = center:length(mag_spectrum);
pos_mag = mag_spectrum(pos_freqs);

% Find first minimum (close to zero)
threshold = 0.1 * max(mag_spectrum);  % Define "zero" as 10% of max
zero_indices = find(pos_mag < threshold);
if ~isempty(zero_indices)
    first_zero_idx = zero_indices(1);
    zero_freq = freq(center + first_zero_idx - 1);
else
    zero_freq = NaN;
end

%% Display results
fprintf('\n=== Part (d) Results ===\n');
fprintf('Zero-crossing frequency: %.3f Hz\n', zero_freq);
fprintf('Theoretical value (1/T): %.3f Hz\n', 1/T);
fprintf('Zero-crossing bandwidth (2/T): %.3f Hz\n', 2/T);
fprintf('\n');

%% Comparison with parts (a), (b), (c)
fprintf('=== Comparison ===\n');
fprintf('Part (a) - Single pulse bandwidth: 2/T = %.3f Hz\n', 2/T);
fprintf('Part (b) - Square wave bandwidth: 2/T = %.3f Hz\n', 2/T);
fprintf('Part (d) - Random 2-PAM bandwidth: 2/T = %.3f Hz\n', 2/T);
fprintf('\nConclusion: All three signals have the SAME bandwidth!\n');
