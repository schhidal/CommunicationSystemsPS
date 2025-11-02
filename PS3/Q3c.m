% --- MATLAB/Octave Script for Problem 3 ---

clear; close all;

% --- 1. Define Frequencies from the Problem ---
f1 = 2000;         % 2 kHz (Baseband bandwidth)
f2 = 1e6;          % 1 MHz (RF Carrier)
f3 = 990e3;        % 990 kHz (Interferer 1)
f4 = 1005e3;       % 1005 kHz (Interferer 2)
f5 = 995e3;        % 995 kHz (Our chosen BPF cutoff)
f6 = 1003e3;       % 1003 kHz (Our chosen BPF cutoff)
f7 = 950e3;        % 950 kHz (Sampling Frequency)
f8 = 50e3;         % 50 kHz (Our chosen IF-to-Baseband mixer freq)
f9 = 10e3;         % 10 kHz (Our chosen final LPF cutoff)
T = 20e-6;         % 20 us (Symbol Period)
M = 19;            % Our calculated downsampling factor
f_symbol = 1/T;    % 50 kHz (Symbol Rate)

% --- 2. Define Simulation Parameters ---
% We need a simulation sampling rate fs_sim that is
% higher than all signals, including interferers and replicas.
% The BPF will filter x3, which has components up to 1005 kHz.
% Let's sample at 2.5 MHz.
fs_sim = 2.5e6;
Ts_sim = 1/fs_sim;
N = 125000;
t_sim = (0:N-1) * Ts_sim;

disp('Parameters set. Simulating at 2.5 MHz.');

% --- 3. Create Baseband Signal x1(t) ---
disp('Generating x1(t)...');
white_noise = randn(1, N);
[b_lpf1, a_lpf1] = butter(6, f1 / (fs_sim/2));
x1 = filter(b_lpf1, a_lpf1, white_noise);
figure;
plotspec(x1, Ts_sim);
title('x_1(t) - Baseband Signal');
subplot(2,1,2); xlim([-f1*3, f1*3]);

% --- 4. Create RF Signal x2(t) ---
disp('Generating x2(t)...');
x2 = x1 .* cos(2*pi*f2*t_sim);

% --- 5. Create "Dirty" Signal x3(t) with Interferers ---
disp('Generating x3(t)...');
interferer1 = 0.5 * cos(2*pi*f3*t_sim);
interferer2 = 0.5 * cos(2*pi*f4*t_sim);
x3 = x2 + interferer1 + interferer2;

% --- 6. Create "Clean" RF Signal x4(t) with BPF ---
disp('Generating x4(t)...');
[b_bpf, a_bpf] = butter(6, [f5, f6] / (fs_sim/2), 'bandpass');
x4 = filter(b_bpf, a_bpf, x3);
figure;
plotspec(x4, Ts_sim);
title('x_4(t) - "Clean" RF Signal (Interferers Removed)');
subplot(2,1,2); xlim([f5-10e3, f6+10e3]);

% --- 7. Sample x4(t) to get x5(nTs) ---
% This is the sub-Nyquist sampling step.
% We simulate this by "downsampling" our simulation signal.
disp('Generating x5(nTs) via sub-sampling...');
M_sub = fs_sim / f7; % M_sub = 2.5e6 / 9.5e5 = 2.63... not integer
% Since it's not integer, we'll "resample" (interpolate and decimate)
% This is a more realistic way to simulate this step.
x5 = resample(x4, f7, fs_sim);
Ts_5 = 1/f7;
figure;
plotspec(x5, Ts_5);
title(['x_5(nTs) - Sub-sampled to ', num2str(f7/1e3), ' kHz (Aliased to +/- 50 kHz)']);
subplot(2,1,2); xlim([-f7/2, f7/2]);

% --- 8. Mix x5 to baseband to get x6(nTs) ---
disp('Generating x6(nTs)...');
t_5 = (0:length(x5)-1) * Ts_5; % New time vector for x5
c_if = cos(2*pi*f8*t_5);
x6 = x5 .* c_if;
figure;
plotspec(x6, Ts_5);
title(['x_6(nTs) - Mixed to Baseband (Replicas at +/- 100 kHz)']);
subplot(2,1,2); xlim([-f7/2, f7/2]);

% --- 9. LPF x6 to get x7(nTs) ---
disp('Generating x7(nTs)...');
[b_lpf2, a_lpf2] = butter(6, f9 / (f7/2));
x7 = filter(b_lpf2, a_lpf2, x6);
figure;
plotspec(x7, Ts_5);
title(['x_7(nTs) - Final Baseband Signal (Sampled at ', num2str(f7/1e3), ' kHz)']);
subplot(2,1,2); xlim([-f9*2, f9*2]);

% --- 10. Downsample x7 to get x8(nT) ---
disp('Generating x8(nT)...');
% We just take every M-th sample
x8 = x7(1:M:end);
Ts_8 = 1/f_symbol;
disp(['Final signal x8 has ', num2str(length(x8)), ' samples.']);
disp('Receiver design complete.');
figure;
plotspec(x8, Ts_8);