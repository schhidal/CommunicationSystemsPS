Ts = 0.001;    % Sampling interval
T = 10;     % Total Time
t = 0:Ts:T-Ts;  % Time vector (length N = T /Ts)


% Signals x1 to x7
x1 = (t >= 0);
x2 = (t >= 0) - (t > 1);
x3 = exp(-t) .* (t >= 0);
x4 = cos (2 * pi * 1 * t) .* (t >= 0);
x5 = (1 - exp(-2 .* t)) .* (exp(-t) .* (t >= 0));
x6 = (cos(2*pi*t) .* (exp(-t) .* (t >= 0)));
x7 = (2*rand(size(t))-1);

figure('Name','Spectrum');
plotspec(x1,Ts);
figure('Name','Spectrum');
plotspec(x2,Ts);
figure('Name','Spectrum');
plotspec(x3,Ts);
figure('Name','Spectrum');
plotspec(x4,Ts);
figure('Name','Spectrum');
plotspec(x5,Ts);
figure('Name','Spectrum');
plotspec(x6,Ts);
figure('Name','Spectrum');
plotspec(x7,Ts);




