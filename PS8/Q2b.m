%% Part (b): Receiver Implementation to Reconstruct Message from m6params.m
% This receiver implements the key components needed to decode the 4-PAM signal

clear all; close all; clc;

%% ===== LOAD TRANSMITTED SIGNAL FROM PART (a) =====
% If you saved results from part (a), load them
% Otherwise, regenerate the signal
if exist('part_a_results.mat', 'file')
    load('part_a_results.mat');
    fprintf('Loaded signal from part (a)\n');
else
    % Regenerate the signal
    fprintf('Regenerating signal...\n');
    m = '0123456789 I wish I were an Oscar Meyer wiener 56789';

    frameParams.userDataLength = 1;
    frameParams.preamble = '';
    frameParams.chanCodingFlag = 0;
    frameParams.bitEncodingFlag = 0;

    rfParams.f_s = 100;
    rfParams.T_t = 1;
    rfParams.T_t_err = 0;
    rfParams.f_if = 20;
    rfParams.f_if_err = 0;
    rfParams.phaseNoiseVariance = 0;
    rfParams.SRRCLength = 4;
    rfParams.SRRCrolloff = 0.3;

    chanParams.c1 = [1 0 0];
    chanParams.c2 = [1 0 0];
    chanParams.randomWalkVariance = 0;
    chanParams.SNR = Inf;
    chanParams.adjacentUser1Power = -Inf;
    chanParams.adjacentUser1f_if = 0;
    chanParams.adjacentUser1Chan = [1 0 0];
    chanParams.adjacentUser2Power = -Inf;
    chanParams.adjacentUser2f_if = 0;
    chanParams.adjacentUser2Chan = [1 0 0];
    chanParams.NBIfreq = 0;
    chanParams.NBIPower = -Inf;

    [r, s] = BigTransmitter(m, frameParams, rfParams, chanParams);
end

fprintf('\n===== RECEIVER IMPLEMENTATION =====\n');
fprintf('Received signal length: %d samples\n', length(r));
fprintf('Expected number of symbols: %d\n', length(s));

%% ===== RECEIVER PARAMETERS =====
M = rfParams.f_s;                         % Oversampling factor (samples per symbol)
f_c = rfParams.f_if;                      % Carrier frequency
T_s = 1/rfParams.f_s;                     % Sampling period
T_symbol = rfParams.T_t;                  % Symbol period

%% ===== STEP 1: DIGITAL DOWNCONVERSION (IF to Baseband) =====
fprintf('\nStep 1: Digital Downconversion...\n');

% Create time vector
t = (0:length(r)-1) * T_s;

% Assume perfect carrier synchronization (ideal case for now)
% In reality, this is the SHOWSTOPPER - we would need carrier recovery here
theta = 0;  % Phase offset (assumed known in ideal case)

% Downconvert using cosine (in-phase) and sine (quadrature)
% For real 4-PAM, we only need the in-phase component
c_cos = 2 * cos(2*pi*f_c*t + theta);

% Mix to baseband
x_i = r(:) .* c_cos(:);  % Ensure column vectors

fprintf('Downconversion complete. Using in-phase component.\n');

%% ===== STEP 2: LOWPASS FILTER (Remove double-frequency components) =====
fprintf('\nStep 2: Lowpass Filtering...\n');

% Design lowpass filter
f_cutoff = 0.1;  % Normalized cutoff frequency
lpf_order = 50;
b_lpf = firpm(lpf_order, [0, f_cutoff, f_cutoff+0.05, 1], [1, 1, 0, 0]);

% Filter the downconverted signal
x_lpf = filter(b_lpf, 1, x_i);

% Account for filter delay
lpf_delay = floor(lpf_order/2);

fprintf('Lowpass filtering complete. Filter delay: %d samples\n', lpf_delay);

%% ===== STEP 3: PULSE-MATCHED FILTER (SRRC Correlator) =====
fprintf('\nStep 3: Matched Filtering with SRRC...\n');

% Generate SRRC pulse for matched filtering
srrc_pulse = srrc(rfParams.SRRCLength, rfParams.SRRCrolloff, M, 1);

% Correlate with SRRC pulse (matched filter)
y = filter(srrc_pulse, 1, x_lpf);
y = y(:);  % Ensure column vector

% Account for matched filter delay
srrc_delay = floor(length(srrc_pulse)/2);
total_delay = lpf_delay + srrc_delay;

fprintf('Matched filtering complete. SRRC delay: %d samples\n', srrc_delay);
fprintf('Total processing delay: %d samples\n', total_delay);
fprintf('Signal length after matched filter: %d samples\n', length(y));

%% ===== STEP 4: TIMING RECOVERY AND DOWNSAMPLING =====
fprintf('\nStep 4: Timing Recovery and Downsampling...\n');

% Compute absolute value for peak detection
abs_y = abs(y);

% Find maximum and set threshold
max_val = max(abs_y);
threshold = 0.5 * max_val;

fprintf('Maximum signal value: %.4f\n', max_val);
fprintf('Threshold for peak detection: %.4f\n', threshold);

% Safe search start (ensure it's within bounds)
search_start = min(total_delay, length(abs_y) - M);
search_start = max(1, search_start);

fprintf('Searching for first peak starting at sample %d\n', search_start);

% Find first peak
if search_start < length(abs_y)
    search_region = abs_y(search_start:end);
    peaks_idx = find(search_region > threshold, 1, 'first');

    if ~isempty(peaks_idx)
        data_start = search_start + peaks_idx - 1;
    else
        % If no peak found, use default
        data_start = total_delay + M;
        fprintf('Warning: No peak found above threshold. Using default start.\n');
    end
else
    data_start = M;
    fprintf('Warning: Search start beyond signal length. Using default.\n');
end

fprintf('Initial data start: %d\n', data_start);

% Fine-tune to land on symbol center (optimize over one symbol period)
% This finds the sampling phase that maximizes power
best_power = 0;
best_offset = 0;

for offset = 0:M-1
    test_start = data_start + offset;
    test_end = min(test_start + 20*M, length(y));

    if test_start <= length(y) && test_end > test_start
        test_indices = test_start:M:test_end;
        test_indices = test_indices(test_indices <= length(y));

        if length(test_indices) >= 5  % Need at least 5 samples
            test_power = sum(abs(y(test_indices)).^2);
            if test_power > best_power
                best_power = test_power;
                best_offset = offset;
            end
        end
    end
end

data_start = data_start + best_offset;

fprintf('Optimized data start: %d (offset +%d)\n', data_start, best_offset);
fprintf('Best power: %.4f\n', best_power);

% Downsample to get one sample per symbol
if data_start <= length(y)
    z = y(data_start:M:end);
else
    error('Data start is beyond signal length. Cannot downsample.');
end

fprintf('Downsampled to %d soft decisions\n', length(z));

%% ===== STEP 5: NORMALIZE (AGC) =====
fprintf('\nStep 5: Automatic Gain Control...\n');

% Calculate normalization factor based on expected symbol power
% For 4-PAM alphabet [-3, -1, 1, 3], expected power is (9+1+1+9)/4 = 5
expected_symbol_power = 5;
actual_power = mean(abs(z).^2);

if actual_power > 0
    receiver_gain = sqrt(expected_symbol_power / actual_power);
else
    receiver_gain = 1;
    fprintf('Warning: Zero signal power detected.\n');
end

% Normalize the soft decisions
z_normalized = z * receiver_gain;

fprintf('Actual signal power: %.4f\n', actual_power);
fprintf('Receiver gain: %.6f\n', receiver_gain);
fprintf('Normalized power: %.4f\n', mean(abs(z_normalized).^2));

%% ===== STEP 6: QUANTIZATION (Decision Device) =====
fprintf('\nStep 6: Quantization...\n');

% Define 4-PAM alphabet
alphabet = [-3, -1, 1, 3];

% Quantize to nearest symbol
m_hat = quantalph(z_normalized, alphabet);

fprintf('Quantized %d symbols\n', length(m_hat));

% Show symbol distribution
for sym = alphabet
    count = sum(m_hat == sym);
    fprintf('  Symbol %+d: %d occurrences (%.1f%%)\n', sym, count, 100*count/length(m_hat));
end

%% ===== STEP 7: DECODING (Symbols to Text) =====
fprintf('\nStep 7: Decoding...\n');

% Make sure we have a multiple of 4 symbols (required for pam2letters)
% Each character requires 4 PAM symbols (2 bits per symbol, 8 bits per char)
remainder = mod(length(m_hat), 4);
if remainder ~= 0
    % Truncate to nearest multiple of 4
    m_hat_truncated = m_hat(1:end-remainder);
    fprintf('Truncating %d symbols to make length divisible by 4\n', remainder);
else
    m_hat_truncated = m_hat;
end

fprintf('Decoding %d symbols into %d characters...\n', ...
    length(m_hat_truncated), length(m_hat_truncated)/4);

% Decode symbols back to text
try
    if frameParams.bitEncodingFlag == 0
        % 8-bit ASCII encoding (letters2pam/pam2letters)
        message_decoded = pam2letters(m_hat_truncated);
    else
        % 7-bit ASCII encoding
        message_decoded = bin2text(pam2letters(m_hat_truncated));
    end
catch ME
    fprintf('Error during decoding: %s\n', ME.message);
    % Alternative decoding if pam2letters fails
    fprintf('Attempting alternative decoding...\n');
    message_decoded = decode_pam_manual(m_hat_truncated);
end

%% ===== RESULTS =====
fprintf('\n===== RECONSTRUCTION RESULTS =====\n');
fprintf('Original message length: %d characters\n', length(m));
fprintf('Decoded message length:  %d characters\n', length(message_decoded));
fprintf('\nOriginal message:\n%s\n', m);
fprintf('\nDecoded message:\n%s\n', message_decoded);

% Calculate error rate
min_len = min(length(m), length(message_decoded));
if min_len > 0
    errors = sum(m(1:min_len) ~= message_decoded(1:min_len));
    error_rate = errors / min_len * 100;

    fprintf('\nCharacter error rate: %.2f%% (%d errors out of %d)\n', ...
        error_rate, errors, min_len);

    % Also calculate symbol error rate
    if length(s) >= length(m_hat)
        s_expected = s(1:length(m_hat));
        symbol_errors = sum(s_expected(:) ~= m_hat(:));
        symbol_error_rate = symbol_errors / length(s_expected) * 100;
        fprintf('Symbol error rate: %.2f%% (%d errors out of %d)\n', ...
            symbol_error_rate, symbol_errors, length(s_expected));
    end
else
    fprintf('Cannot calculate error rate - no decoded message\n');
end

%% ===== VISUALIZATION =====
figure('Position', [100 100 1400 900]);

% Plot 1: Received signal
subplot(4,2,1);
plot_len = min(5000, length(r));
plot(t(1:plot_len), r(1:plot_len), 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal r(t)');
grid on;

% Plot 2: Downconverted signal (before LPF)
subplot(4,2,2);
plot_len = min(5000, length(x_i));
plot((0:plot_len-1)*T_s, x_i(1:plot_len), 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Downconverted Signal (before LPF)');
grid on;

% Plot 3: After lowpass filter
subplot(4,2,3);
plot_len = min(5000, length(x_lpf));
plot((0:plot_len-1)*T_s, x_lpf(1:plot_len), 'g');
xlabel('Time (s)');
ylabel('Amplitude');
title('After Lowpass Filter');
grid on;

% Plot 4: After matched filter
subplot(4,2,4);
plot(y, 'k', 'LineWidth', 0.5);
xlabel('Sample index');
ylabel('Amplitude');
title('After Matched Filter');
grid on;
hold on;
% Mark sampling points
sample_indices = data_start:M:min(data_start+50*M, length(y));
if ~isempty(sample_indices)
    plot(sample_indices, y(sample_indices), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
end
legend('Filtered signal', 'Sampling points');

% Plot 5: Eye diagram
subplot(4,2,5);
eye_length = 2*M;  % Two symbols wide
num_traces = min(100, floor((length(y)-data_start)/M) - 2);
hold on;
for i = 1:num_traces
    start_idx = data_start + (i-1)*M;
    if start_idx + eye_length - 1 <= length(y)
        plot(0:eye_length-1, y(start_idx:start_idx+eye_length-1), 'b', 'LineWidth', 0.5);
    end
end
xlabel('Sample index');
ylabel('Amplitude');
title('Eye Diagram');
grid on;

% Plot 6: Constellation (soft decisions)
subplot(4,2,6);
plot(real(z_normalized), imag(z_normalized), 'b.', 'MarkerSize', 8);
hold on;
plot(alphabet, zeros(size(alphabet)), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation Diagram (Soft Decisions)');
grid on;
legend('Soft decisions', 'Alphabet');
axis equal;
ylim([-1 1]);
xlim([-4 4]);

% Plot 7: Soft decisions over time
subplot(4,2,7);
plot_len = min(200, length(z_normalized));
plot(1:plot_len, z_normalized(1:plot_len), 'b.-', 'MarkerSize', 6);
hold on;
plot([1 plot_len], [3 3], 'r--', 'LineWidth', 1);
plot([1 plot_len], [1 1], 'r--', 'LineWidth', 1);
plot([1 plot_len], [-1 -1], 'r--', 'LineWidth', 1);
plot([1 plot_len], [-3 -3], 'r--', 'LineWidth', 1);
plot([1 plot_len], [0 0], 'k:', 'LineWidth', 0.5);
xlabel('Symbol index');
ylabel('Amplitude');
title('Soft Decisions Over Time');
grid on;
legend('Soft decisions', 'Decision boundaries');

% Plot 8: Hard decisions (reconstructed symbols)
subplot(4,2,8);
num_to_show = min(100, length(m_hat));
stem(1:num_to_show, m_hat(1:num_to_show), 'filled', 'MarkerSize', 4);
xlabel('Symbol index');
ylabel('Symbol value');
title(sprintf('Hard Decisions (first %d symbols)', num_to_show));
grid on;
ylim([-4 4]);
yticks([-3 -1 1 3]);

%% ===== SAVE RESULTS =====
save('part_b_results.mat', 'message_decoded', 'm_hat', 'z_normalized');
fprintf('\nResults saved to part_b_results.mat\n');

%% ===== HELPER FUNCTION: MANUAL PAM DECODER =====
function text_out = decode_pam_manual(symbols)
% Manual decoder for 4-PAM to ASCII
% Each 4 symbols (2 bits each) = 8 bits = 1 character

% Ensure length is multiple of 4
n_symbols = length(symbols);
n_chars = floor(n_symbols / 4);
symbols = symbols(1:n_chars*4);

% Map symbols to bits: -3→00, -1→01, 1→10, 3→11
bits = zeros(1, n_symbols * 2);
for i = 1:n_symbols
    if symbols(i) == -3
        bits(2*i-1:2*i) = [0 0];
    elseif symbols(i) == -1
        bits(2*i-1:2*i) = [0 1];
    elseif symbols(i) == 1
        bits(2*i-1:2*i) = [1 0];
    else  % symbols(i) == 3
        bits(2*i-1:2*i) = [1 1];
    end
end

% Convert bits to characters
text_out = char(zeros(1, n_chars));
for i = 1:n_chars
    byte_bits = bits((i-1)*8+1:i*8);
    ascii_val = sum(byte_bits .* (2.^(7:-1:0)));
    text_out(i) = char(ascii_val);
end
end
