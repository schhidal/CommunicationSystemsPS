% BigIdeal_PartC.m
% Solves Part (c): Carrier Phase Offset and AGC Limit

clear; clc; close all;

% --- 1. TRANSMITTER SETUP (Standard) ---
m = '01234 I wish I were an Oscar Meyer wiener 56789';
frameParams.userDataLength = length(m); 
frameParams.preamble = '';
frameParams.chanCodingFlag = 0;
frameParams.bitEncodingFlag = 0;

% Ideal Channel
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

[r, s_trans] = BigTransmitter(m, frameParams, rfParams, chanParams);


% --- 2. RECEIVER WITH PHASE OFFSET ---
M = 100;
fc = 20;
srrc_length = 4;
data_start = 5;

% *** INTRODUCE PHASE OFFSET HERE ***
% Try changing this: 0, 0.5, 1.0, 1.4, 1.5 (approx pi/2)
phi_offset = 1.0; % Radians (approx 57 degrees)

t = (0:length(r)-1)/rfParams.f_s;
% Receiver uses cos(2*pi*fc*t + phi) where phi is the error
c2 = cos(2*pi*fc*t' + phi_offset); 

x2 = r .* c2; 

% Match Filter
filt_pulse = srrc(srrc_length, rfParams.SRRCrolloff, M);
y = filter(filt_pulse, 1, x2);

% Downsample
z_raw = y(srrc_length*M : M : end);
z_raw = z_raw(data_start : data_start + length(s_trans) - 1);


% --- 3. AGC LOOP (From Part b) ---
g = zeros(size(z_raw)); 
g(1) = 1;
mu = 0.001;
ds = 5;
z_agc = zeros(size(z_raw));

for k = 1 : length(z_raw) - 1
    z_agc(k) = g(k) * z_raw(k);
    err = z_agc(k)^2 - ds;
    g(k+1) = g(k) - mu * err; 
end
z_agc(end) = g(end) * z_raw(end);


% --- 4. ANALYSIS ---
fprintf('Phase Offset: %.2f rad (%.1f deg)\n', phi_offset, phi_offset*180/pi);
fprintf('Theoretical Attenuation cos(phi): %.4f\n', cos(phi_offset));
fprintf('Final AGC Gain: %.4f\n', g(end));

% Check decoding
mprime = quantalph(z_agc, [-3, -1, 1, 3])';
errors = sum(abs(sign(mprime - s_trans)))/2;

if errors == 0
    fprintf('Status: SUCCESS (Message recovered)\n');
else
    fprintf('Status: FAILED (%d errors)\n', errors);
end


% --- 5. PLOTS ---
figure;
subplot(2,1,1);
plot(g);
title(['AGC Gain (Phase Offset \phi = ' num2str(phi_offset) ')']);
xlabel('Symbol Index'); ylabel('Gain');

subplot(2,1,2);
plot(z_agc, '.');
hold on;
% Octave compatible lines
levels = [-3 -1 1 3];
for val = levels
    plot([1 length(z_agc)], [val val], 'k:');
end
title('Soft Decisions');
ylim([-5 5]);

% decode decision device output to text string
reconstructed_message=pam2letters(mprime)   % reconstruct message