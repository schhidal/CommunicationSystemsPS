% --- Part (c): Improve SNR with FIR Filter (Corrected Plotting) ---

clear; clc; close all;

% --- Function for Average Power ---
function p = pow(vec)
    p = sum(vec .* conj(vec)) / length(vec); % Use conj for complex compatibility
end

% --- Setup ---
N_samples = 200; % Number of samples for signals
noise_gain = 0.5; % Adjust to control noise level

% Define the impulse response h[n] = delta[n-1] - delta[n-2]
h = [0, 1, -1];
n_h = 0:(length(h)-1);

% Define the unit step input u[n]
n_u = 0:(N_samples-1);
u_n = ones(1, N_samples);

% Define white noise w[n]
w_n = randn(1, N_samples);

% --- Create Noisy Input ---
x_n = u_n + noise_gain * w_n;

% --- Output of Original System h[n] ---
y_noisy = conv(x_n, h);
% Output length and index
N_y = length(y_noisy);
n_y = 0:(N_y-1);

% Separate signal and noise components *before* LPF
s_out = conv(u_n, h);
n_out = conv(noise_gain * w_n, h);
% Need to pad s_out and n_out to match y_noisy length
s_out_padded = [s_out, zeros(1, N_y - length(s_out))]; % Use padded versions for calculations/plotting
n_out_padded = [n_out, zeros(1, N_y - length(n_out))];

% --- Calculate SNR Before LPF ---
% Use trimmed versions for power calculation if desired, or padded ok if long enough
SNR_before = pow(s_out_padded) / pow(n_out_padded);
fprintf('SNR before LPF: %.2f (%.2f dB)\n', SNR_before, 10*log10(SNR_before));

% --- Design the FIR LPF ---
N_taps_lpf = 51; % Number of taps for the LPF
filter_order_lpf = N_taps_lpf - 1;
cutoff_freq_norm = 0.25; % Normalized cutoff frequency (e.g., pi/2 -> 0.25*2pi)
trans_width_norm = 0.1;  % Normalized transition width

f_lpf_norm = [0, cutoff_freq_norm, cutoff_freq_norm + trans_width_norm, 1];
a_lpf = [1, 1, 0, 0];

b_lpf = firpm(filter_order_lpf, f_lpf_norm, a_lpf);

% --- Apply LPF ---
y_final = filter(b_lpf, 1, y_noisy);
N_final = length(y_final); % y_final will have same length as y_noisy
n_final = 0:(N_final-1); % Correct index for y_final

% Apply LPF to separated components to calculate output SNR
s_final = filter(b_lpf, 1, s_out_padded); % Filter the padded signal component
n_final = filter(b_lpf, 1, n_out_padded); % Filter the padded noise component

% --- Calculate SNR After LPF ---
SNR_after = pow(s_final) / pow(n_final);
fprintf('SNR after LPF:  %.2f (%.2f dB)\n', SNR_after, 10*log10(SNR_after));

% --- Plotting Results ---
figure('Name', 'SNR Improvement with LPF');

subplot(3, 1, 1);
stem(n_y, y_noisy, 'filled', 'MarkerSize', 4); hold on;
signal_peak_n_before = 1; % Expected peak index n=1 for delta[n-1]
plot(signal_peak_n_before, 1, 'rx', 'LineWidth', 2, 'MarkerSize', 10); % Use plot for single marker
hold off;
title('Output of h[n] (Signal \delta[n-1] + Filtered Noise)');
xlabel('n'); ylabel('Amplitude'); xlim([-2 N_samples+5]); grid on;
legend('y_{noisy}[n]', 'Original Signal Pulse Location');

subplot(3, 1, 2);
stem(n_final, y_final, 'filled', 'MarkerSize', 4); hold on; % Plot the final filtered output
% Calculate expected delay of LPF
lpf_delay = filter_order_lpf / 2;
signal_peak_n_after = 1 + lpf_delay; % Expected n index of the peak after LPF delay

% Find amplitude around the expected peak time for plotting marker height
peak_amplitude_approx = max(s_final); % Use max value of filtered signal component as Y value

% Check if the expected peak is within the plotted range before plotting
if signal_peak_n_after >= min(n_final) && signal_peak_n_after <= max(n_final)
    plot(signal_peak_n_after, peak_amplitude_approx, 'rx', 'LineWidth', 2, 'MarkerSize', 10); % Mark expected peak
end
hold off;
title('Output After Applying LPF');
xlabel('n'); ylabel('Amplitude'); xlim([-2 N_samples+5+lpf_delay]); grid on;
legend('y_{final}[n]', 'Expected Signal Pulse Location');

% Optional: Plot Filtered Noise Spectrum (ensure plotspec.m is available)
% figure; plotspec(n_out_padded,1); title('Spectrum of Noise Before LPF (High-Pass Char.)');
% figure; plotspec(n_final,1); title('Spectrum of Noise After LPF');