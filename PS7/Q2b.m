% BigIdeal_PartB.m
% Solves Part (b): Adding Automatic Gain Control (AGC)

clear; clc; close all;
addpath('../')

% --- 1. TRANSMITTER & CHANNEL SETUP (Same as Part a) ---
m = '01234 I wish I were an Oscar Meyer wiener 56789';
frameParams.userDataLength = length(m);
frameParams.preamble = '';
frameParams.chanCodingFlag = 0;
frameParams.bitEncodingFlag = 0;

% Ideal Channel (High SNR)
chanParams.c1 = [1 0 0];
chanParams.c2 = [1 0 0];
chanParams.randomWalkVariance = 0;
chanParams.SNR = 100;
chanParams.adjacentUser1Power = -Inf;
chanParams.adjacentUser1f_if = 0;
chanParams.adjacentUser1Chan = [1 0 0];
chanParams.adjacentUser2Power = -Inf;
chanParams.adjacentUser2f_if = 0;
chanParams.adjacentUser2Chan = [1 0 0];
chanParams.NBIfreq = 0;
chanParams.NBIPower = -Inf;

% RF Parameters
rfParams.f_if_err = 0;
rfParams.T_t_err = 0;
rfParams.phaseNoiseVariance = 0;
rfParams.SRRCLength = 4;
rfParams.SRRCrolloff = 0.3;
rfParams.f_s = 100;
rfParams.T_t = 1;
rfParams.f_if = 20;

% Generate Signal
[r, s_trans] = BigTransmitter(m, frameParams, rfParams, chanParams);

% --- 2. RECEIVER FRONT END ---
M = 100;
fc = 20;
srrc_length = 4;
data_start = 5;

% Demodulate
t = (0:length(r)-1)/rfParams.f_s;
x2 = r .* cos(2*pi*fc*t');

% Match Filter
filt_pulse = srrc(srrc_length, rfParams.SRRCrolloff, M);
y = filter(filt_pulse, 1, x2);

% Downsample to get raw (unscaled) soft decisions
z_raw = y(srrc_length*M : M : end);
z_raw = z_raw(data_start : data_start + length(s_trans) - 1);


% --- 3. PART (a) CALCULATION (For Comparison) ---
P_z = mean(z_raw.^2);
P_target = 5; % Average power of {-3, -1, 1, 3}
ideal_fixed_gain = sqrt(P_target / P_z);
fprintf('Part (a) Empirical Fixed Gain: %.4f\n', ideal_fixed_gain);


% --- 4. PART (b) AUTOMATIC GAIN CONTROL ---
% Initialize AGC variables
g = zeros(size(z_raw));
g(1) = 1;        % Start with unity gain (or a guess)
mu = 0.001;      % Step size (controls convergence speed)
ds = 5;          % Desired Signal Power (Target)
z_agc = zeros(size(z_raw));

for k = 1 : length(z_raw) - 1
    % 1. Apply current gain
    z_agc(k) = g(k) * z_raw(k);

    % 2. Measure error (Power difference)
    %    output_power - target_power
    err = z_agc(k)^2 - ds;

    % 3. Update gain (Gradient Descent)
    %    If power is too high (err > 0), decrease gain.
    g(k+1) = g(k) - mu * err;
end
% Apply final gain to last symbol
z_agc(end) = g(end) * z_raw(end);

fprintf('Part (b) Final AGC Gain:       %.4f\n', g(end));


% --- 5. PLOTS & VERIFICATION ---
figure;
subplot(2,1,1);
plot(g);
hold on;
% Octave-compatible replacement for yline
plot([1 length(g)], [ideal_fixed_gain ideal_fixed_gain], 'r--');
title('AGC Gain Trajectory');
xlabel('Symbol Index'); ylabel('Gain Value');
legend('Adaptive Gain', 'Fixed Target');

subplot(2,1,2);
plot(z_agc, '.');
hold on;
% Octave-compatible replacement for yline([-3 -1 1 3])
levels = [-3 -1 1 3];
for val = levels
    plot([1 length(z_agc)], [val val], 'k:');
end
title('Soft Decisions after AGC');
ylim([-5 5]);

% Decode
mprime = quantalph(z_agc, [-3, -1, 1, 3])';
str = pam2letters(mprime);
fprintf('Recovered Message: "%s"\n', str);

