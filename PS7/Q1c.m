% Modifications for Part (a)

clear; clc; close all;
addpath('../')

% --- TRANSMITTER (Standard idsys.m code) ---
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2pam(str); N=length(m);
M=100;
mup=zeros(1,N*M);
mup(1:M:N*M)=m;
p=hamming(M);
x=filter(p,1,mup);
% figure(1), plotspec(x,1/M) % Commented out to focus on results
t=1/M:1/M:length(x)/M;
fc=20;
c=cos(2*pi*fc*t);
r=c.*x;

%CHANNEL - ADD GAUSSIAN NOISE
% Define noise gain (controls noise power)
noise_gain = 0.6;  % Try 0 (no noise), 0.6 (mild), 2 (harsh)
% Generate white Gaussian noise with same size as received signal
noise = noise_gain * randn(size(r));

% Add noise to received signal
r = r + noise;

% Optional: Calculate and display SNR
signal_power = mean(r.^2);
noise_power = mean(noise.^2);
SNR_dB = 10*log10(signal_power/noise_power);
fprintf('Signal-to-Noise Ratio: %.2f dB\n', SNR_dB);

figure;
plotspec(r,1/M);
subplot(2,1,2);
title('Spectrum after adding Noise');

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
x2=r.*c2;                      % demod received signal
fl=50; fbe=[0 0.1 0.2 1];      % LPF parameters
damps=[1 1 0 0 ];
b=firpm(fl,fbe,damps);         % create LPF impulse response
x3=2*filter(b,1,x2);           % LPF and scale signal
% extract upsampled pulses using correlation implemented
% as a convolving filter; filter with pulse and normalize
y=filter(fliplr(p)/(pow(p)*M),1,x3);
% set delay to first symbol-sample and increment by M
z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
% Standard decoding steps (will likely fail)
mprime=quantalph(z,[-3,-1,1,3])';
cvar=(mprime-z)*(mprime-z)'/length(mprime);  % cluster variance
pererr=100*sum(abs(sign(mprime-m(1:length(mprime)))))/length(mprime);

reconstructed_message=pam2letters(mprime)

% Soft decisions show noise spreading
figure;
plot(1:length(z), z, '.')
hold on
plot(1:length(z), mprime, 'ro', 'MarkerSize', 8)
title(['Soft Decisions with Noise, Gain = ', num2str(noise_gain)])
xlabel('Symbol Index')
ylabel('Amplitude')
legend('Soft Decisions', 'Hard Decisions')
grid on

% Shows distribution around symbol values
figure;
hist(z, 50)
title(['Distribution of Soft Decisions, Noise Gain = ', num2str(noise_gain)])
xlabel('Amplitude')
ylabel('Count')
grid on
hold on

