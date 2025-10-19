Ts = 0.1;
T = 10;
t = 0:Ts:T-Ts;

f_0 = [1, 2, 6, 12];   % fundamental frequencies Q2b
phases =  [0, pi/4, pi/2, pi];

 
for i = 1:length(f_0)
    f_new = f_0(i);
    for j = 1:length(phases);

    phase = phases(i);
    x = cos(2*pi*f_new*t+phase);

figure('Name', 'Spectrum of Sinusoids');
plotspec(x, Ts);

    end
end 