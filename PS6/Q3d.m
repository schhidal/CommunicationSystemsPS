% ## Partb

%Short 4-PAM message 
m= letters2pam('Hello');    % Convert to 4-PAM symbols {-3, -1, 1, +3}

% Parameters 
M= 100;             % Oversampling factor (samples per symbol)
T=1;                %# Symbol period 


% Generate pulse shape: Hamming blip 
ps= hamming(M);


% Oversample the sysmbols 
mup = zeros(1, length(m)*M); 
mup(1:M:end)=m;         % Place symbols at M intervals

% Pulse shaping (convolve) 
x= filter(ps, 1,mup);   %# Transmitted pulse shaped signal 

% Inspect signal and spectrum 
figure;
plotspec(x, 1/M);
title('Pulse-Shaped 4-PAM signal with Hamming Blip');



% Ajustable noise power parameter
noise_power= 0.2;

% Generate white noise 
n = sqrt (noise_power) * randn(size(x));

% Add noise to signal 
x_noisy = x + n;

% Plot noisy signal and its spectrum 
figure;
plotspec(x_noisy,1/M);
title('Noisy Pulse-Shaped Signal');


% part c
% Correlate with the Hamming blip 
% xcorr performs correlation 
y = xcorr(x_noisy, ps);

% Normalize by pulse energy 

pulse_energy = sum(ps.^2);
y_normalized = y/ pulse_energy;

% plot the correlation result 

figure;
plotspec(y_normalized, 1/M);
subplot(2,1,2);
title('Spectrum After Correlation')



% Part d 

%% Part (d): Correlation vs Convolution

% Method 1: Using xcorr (from part c)
y_xcorr = xcorr(x_noisy, ps);
y_xcorr = y_xcorr / sum(ps.^2);

% Method 2: Using conv
y_conv = conv(x_noisy, fliplr(ps));
y_conv = y_conv / sum(ps.^2);

% Method 3: Using filter (matched filter)
y_filt = filter(fliplr(ps), 1, x_noisy);
y_filt = y_filt / sum(ps.^2);

% Compare
figure;
subplot(3,1,1); plot(y_xcorr); title('Correlation (xcorr)');
subplot(3,1,2); plot(y_conv); title('Convolution (conv)');
subplot(3,1,3); plot(y_filt); title('Filter (filter)');
