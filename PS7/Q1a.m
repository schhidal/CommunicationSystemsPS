% --- idsys_freq_offset.m ---
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

% --- RECEIVER (Modified for Frequency Offset) ---

% 1. Define the frequency offset gamma
gamma = 0.01; % Try 0.01, 0.05, 0.1 Hz

% 2. Create the receiver's unsynchronized oscillator
% The receiver thinks the freq is fc, but effectively it's creating
% a beat frequency because of the mismatch.
c2 = cos(2*pi*(fc + gamma)*t);

% 3. Demodulate
x2 = r .* c2;

% 4. Low Pass Filter (Standard)
fl=50; fbe=[0 0.1 0.2 1];
damps=[1 1 0 0 ];
b=firpm(fl,fbe,damps);
x3=2*filter(b,1,x2);

% 5. Correlate and Downsample (Standard)
y=filter(fliplr(p)/(pow(p)*M),1,x3);
z=y(0.5*fl+M:M:N*M);

% --- VISUALIZATION ---
% This is the plot best suited to characterize the impairment.
figure;
plot(1:length(z), z, '.');
title(['Soft Decisions with Frequency Offset \gamma = ', num2str(gamma), ' Hz']);
xlabel('Symbol Index'); ylabel('Amplitude');
grid on;

% Standard decoding steps (will likely fail)
mprime=quantalph(z,[-3,-1,1,3])';
pererr=100*sum(abs(sign(mprime-m(1:length(mprime)))))/length(mprime);
disp(['Symbol Error Rate: ', num2str(pererr), '%']);

reconstructed_message=pam2letters(mprime)