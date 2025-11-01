% --MATLAB/Octave Script for part(f): IF to BAseband --
fc = 1000;    % Carrier frequency in Hz
f_if = 200;    % Target INtermediate frequency  in Hz
f_baseband_max = 50;    % Max frequency of our Baseband sinc signal in Hz

% -- Define Simulation Parameters --
fs = 5000; % Samling frequency in Hz
Ts = 1/fs; % Sampling Time interval
t = 0:Ts:4; % Time vector from 0 to 4 seconds

disp('Parameters set. Sampling at 5000 Hz.');

% -- Create the original signal RF from y(t) part 1c
disp('Generating RF signal y(t)...');
x_baseband = sinc(100*(t-2));  % Baseband signal
c_rf = cos (2*pi*fc*t); % 1000 Hz Carrier
y_rf = x_baseband .* c_rf;  %  Modulated RF Signal

% -- Create the original signal RF from y(t) part 1c
disp('Generating RF signal y(t)...');
x_baseband = sinc(100*(t-2));  % Baseband signal
c_rf = cos (2*pi*fc*t); % 1000 Hz Carrier
y_rf = x_baseband .* c_rf;  %  Modulated RF Signal

% -- Perform Low side injection Downconversion --
disp('Performing Low-side injection...');
f_lo_low = fc - f_if; % 1000-800 =200 Hz
c_lo = cos(2*pi*f_lo_low*t);  % 800 Hz local oscillator
w_lo = y_rf .* c_lo;  % The new downconverted signal

% --Step 1 : Filter w_lo(t) from part(d) ---
% We must first remove the replicas at +/- 1800 Hz.
% LPF cutoff must be >200+50 Hz and 1800-50 Hz.
% Let's set cutoff to 300 Hz.\
disp('Filtering to create clean IF signal at +|- 200 Hz.');
f_cutoff_1 =300;
% 'butter' designs a Butterworth filter. [b,a] are the Filter coefficients.
[b1, a1] = butter(6, f_cutoff_1 / (fs/2));
y_if = filter(b1, a1, w_lo);


% Plot the spectrum of the IF signal
figure;
plotspec(y_if, Ts);
subplot(2,1,2);
title('1(f) Clean IF Signal (Replicas at +/- 1800 Hz removed)');
xlim([-1000, 1000]); % Zoom in to the Center

% --Step 2 : Mix the IF signal down to baseband (0 Hz)--
disp('Mixing IF signal down to baseband...');
c_baseband = cos(2*pi*f_if*t); % 200 Hz LO
y_mixed = y_if .* c_baseband; % This signal has boxes at 0 Hz and +/- 400 Hz

% Plot the spectrum of the mixed signal
figure;
plotspec(y_mixed, Ts);
subplot(2,1,2);
title('(f) Step 2: Mixed to Baseband (Boxes at 0 Hz and +/- 400 Hz)');
xlim([-1000, 1000]);

disp('Filtering to get final reconstructed signal....');
f_cutoff_2 = 100;
[b2, a2] = butter(6, f_cutoff_2 / (fs/2));
x_reconstructed = filter(b2, a2, y_mixed);

% Plot the spectrum of the final reconstructed signal
figure;
plotspec(x_reconstructed, Ts);
subplot(2,1,2);
title('(f) Step 3: Reconstructed Baseband Signal (Box at 0 Hz)');
xlim([-200, 200]);  % Zoom in to Baseband



