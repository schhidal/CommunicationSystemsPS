addpath('../Tests');

% pulrecsig_with_preprocessing.m
% Combines pulrecsig.m and pllpreprocess.m

%% Part 1: Generate Signal (from pulrecsig.m)
N=10000; M=20; Ts=.0001;             
time=Ts*N*M; t=Ts:Ts:time;           
m=pam(N,4,5);                        
mup=zeros(1,N*M); mup(1:M:N*M)=m;    
ps=hamming(M);                       
s=filter(ps,1,mup);                  
fc=1000; phoff=-1.0;                 % carrier freq. and phase
c=cos(2*pi*fc*t+phoff);              
rsc=s.*c;                            % SUPPRESSED carrier signal

%% Part 2: Preprocessing (from pllpreprocess.m)
r = rsc;                             % suppressed carrier signal
q = r.^2;                            % square the signal

% Design BPF centered at 2*fc
fl = 500;                            % filter length
ff = [0 .38 .39 .41 .42 1];         % freq bands (normalized to 0-1)
fa = [0 0 1 1 0 0];                 % desired amplitudes
h = firpm(fl, ff, fa);              % design BPF

rp = filter(h, 1, q);               % apply BPF -> preprocessed signal

%% Part 3: Estimate Frequency and Phase
fftrBPF = fft(rp);                  % spectrum of preprocessed signal
[m, imax] = max(abs(fftrBPF(1:end/2)));  % find peak
ssf = (0:length(rp)-1)/(Ts*length(rp));  % frequency vector
freqS = ssf(imax);                  % estimated frequency
phasep = angle(fftrBPF(imax));      % estimated phase (includes BPF shift)

% Calculate BPF phase shift at the peak frequency
[IR, ffreqz] = freqz(h, 1, length(rp), 1/Ts);
[~, im] = min(abs(ffreqz - freqS));
phaseBPF = angle(IR(im));           % phase shift of BPF at peak freq

% Correct for BPF phase shift
phaseS = mod(phasep - phaseBPF, pi); % estimated carrier phase

fprintf('\n=== RESULTS ===\n');
fprintf('Original: fc = %d Hz, phoff = %.3f rad (%.2f deg)\n', ...
        fc, phoff, phoff*180/pi);
fprintf('Preprocessed signal at 2*fc:\n');
fprintf('  freqS = %.2f Hz (should be ≈ %d Hz)\n', freqS, 2*fc);
fprintf('  phasep = %.3f rad (raw, includes BPF shift)\n', phasep);
fprintf('  phaseBPF = %.3f rad (BPF phase shift)\n', phaseBPF);
fprintf('Carrier estimate (after correction):\n');
fprintf('  Carrier freq = %.2f Hz (freqS/2)\n', freqS/2);
fprintf('  Carrier phase = %.3f rad (%.2f deg)\n', phaseS, phaseS*180/pi);
fprintf('  Expected phase = %.3f rad (%.2f deg)\n', phoff, phoff*180/pi);
