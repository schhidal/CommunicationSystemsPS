% BigIdeal_PartE_F.m
% Solves Part (e) and (f): Changing Sampling Frequency

clear; clc; close all;

% --- 1. COMMON SETUP ---
m = '01234 I wish I were an Oscar Meyer wiener 56789';
frameParams.userDataLength = length(m); 
frameParams.preamble = '';
frameParams.chanCodingFlag = 0;
frameParams.bitEncodingFlag = 0;

chanParams.c1 = [1 0 0];
chanParams.c2 = [1 0 0];
chanParams.randomWalkVariance = 0;
chanParams.SNR = 100; % Clean channel
chanParams.adjacentUser1Power = -Inf;
chanParams.adjacentUser1f_if = 0;
chanParams.adjacentUser1Chan = [1 0 0];
chanParams.adjacentUser2Power = -Inf;
chanParams.adjacentUser2f_if = 0;
chanParams.adjacentUser2Chan = [1 0 0];
chanParams.NBIfreq = 0;
chanParams.NBIPower = -Inf;

rfParams.f_if_err = 0;
rfParams.T_t_err = 0;
rfParams.phaseNoiseVariance = 0;
rfParams.SRRCLength = 4;
rfParams.SRRCrolloff = 0.3;
rfParams.T_t = 1;
rfParams.f_if = 20; % fc = 20 Hz


% --- 2. TEST CASE LOOP ---
% We will test two sampling frequencies
fs_list = [50, 15]; 

for test_idx = 1:length(fs_list)
    current_fs = fs_list(test_idx);
    fprintf('\n========================================\n');
    fprintf('TESTING SAMPLING FREQUENCY fs = %d Hz\n', current_fs);
    
    rfParams.f_s = current_fs;
    
    % --- Transmit ---
    [r, s_trans] = BigTransmitter(m, frameParams, rfParams, chanParams);
    
    % --- Adaptive Receiver Setup ---
    % 1. Determine Oversampling Factor M
    %    M = fs * T_symbol
    M = current_fs * rfParams.T_t;
    fprintf('  -> Calculated Oversampling Factor M = %d\n', M);
    
    % 2. Check Nyquist (Passband)
    %    Max freq component approx fc + (1+beta)/(2T)
    B_one_sided = (1 + rfParams.SRRCrolloff) / (2 * rfParams.T_t);
    f_max = rfParams.f_if + B_one_sided;
    Nyquist_rate = 2 * f_max;
    
    fprintf('  -> Signal Max Freq: %.2f Hz\n', f_max);
    fprintf('  -> Required Nyquist Rate: %.2f Hz\n', Nyquist_rate);
    
    if current_fs < Nyquist_rate
        fprintf('  -> WARNING: Sub-Nyquist Sampling! Aliasing expected.\n');
    else
        fprintf('  -> Sampling is sufficient.\n');
    end
    
    
    % --- Receiver Processing ---
    fc = rfParams.f_if;
    srrc_length = rfParams.SRRCLength;
    
    % Demodulate
    t = (0:length(r)-1)/current_fs;
    x2 = r .* cos(2*pi*fc*t');
    
    % Match Filter (Pulse Shape)
    % IMPORTANT: The pulse shape must match the new M
    filt_pulse = srrc(srrc_length, rfParams.SRRCrolloff, M);
    y = filter(filt_pulse, 1, x2);
    
    % Downsample
    % Use correlation to find exact start index (Robust Sync)
    known_pattern = s_trans(1:10);
    pat_up = zeros(1, 10*M);
    pat_up(1:M:end) = known_pattern;
    ref_sig = filter(filt_pulse, 1, pat_up);
    
    c = xcorr(y(1:min(length(y), 100*M)), ref_sig);
    [~, peak_idx] = max(abs(c));
    
    % The peak of cross-correlation corresponds to alignment.
    % If lag is 0, signals align at index 1.
    % xcorr returns lags from -(N-1) to (N-1).
    % We need to map this back to array index.
    % Let's use a simpler approach: Find max of convolution-like output
    % or just assume the delay logic from Chapter 8 holds:
    % Delay = filter_delay (TX) + filter_delay (RX) = 2 * group_delay
    % Group delay of linear phase FIR is (N-1)/2.
    % Here length is approx 2*L*M. So delay is approx L*M.
    % Total delay approx 2*L*M.
    
    % Let's try a search around expected delay
    expected_delay = 2 * srrc_length * M; 
    % Search for max energy point in a window around expected delay
    search_window = expected_delay - M : expected_delay + M;
    % Ensure indices are valid
    search_window = search_window(search_window > 0 & search_window < length(y));
    
    [~, best_offset] = max(abs(y(search_window)));
    start_idx = search_window(best_offset);
    
    % Actually, for baseband PAM with SRRC, the peak value is the symbol.
    % So simply finding the max amplitude in the first symbol period works.
    % But we have a delay. 
    
    % Let's trust the calculated index based on filter length:
    % Start sampling exactly at the peak of the pulse.
    % Filter length is 2*srrc_length*M + 1. Group delay is srrc_length*M.
    % Total delay = TX delay + RX delay = srrc_length*M + srrc_length*M = 2*srrc_length*M.
    % Add 1 for 1-based indexing.
    start_idx = 2 * srrc_length * M + 1;

    z_raw = y(start_idx : M : end);
    
    % Ensure length match
    len_to_take = min(length(z_raw), length(s_trans));
    z_raw = z_raw(1:len_to_take);
    s_ref = s_trans(1:len_to_take);
    
    % AGC (Fast)
    g = 1; mu = 0.01; ds = 5;
    z_agc = zeros(size(z_raw));
    for k = 1:length(z_raw)
        z_out = g * z_raw(k);
        z_agc(k) = z_out;
        g = g - mu * (z_out^2 - ds);
    end
    
    % Errors
    mprime = quantalph(z_agc, [-3, -1, 1, 3])';
    errs = sum(abs(sign(mprime - s_ref)))/2;
    
    fprintf('  -> Symbol Errors: %d\n', errs);
    
    % Plot
    figure;
    plot(z_agc, '.');
    hold on;
    % Octave-compatible lines
    for val = [-3 -1 1 3]
        plot([1 length(z_agc)], [val val], 'k:');
    end
    title(sprintf('fs = %d Hz (M=%d)', current_fs, M));
    ylim([-5 5]);
    
    if errs > length(s_trans)/4
         fprintf('  -> STATUS: FAILED (High Error Rate)\n');
         if current_fs == 15
             fprintf('     Explanation: Aliasing caused spectrum overlap.\n');
         end
    else
         fprintf('  -> STATUS: SUCCESS\n');
    end
end