%% After running idsys.m main code...

% EYE DIAGRAM (after pulse correlation)
figure(1);
% Determine delay and number of 4-symbol groups
delay = 0.5*fl + M;  % LPF delay + pulse shape delay
u = floor((length(y)-delay)/(4*M));

% Reshape and plot overlays
plot(reshape(y(delay+1:u*4*M+delay), 4*M, u));
xlabel('Sample within 4-symbol period');
ylabel('Amplitude');
title('Eye Diagram - 4-PAM Signal');
grid on;

% CONSTELLATION HISTORY (after downsampling)
figure(2);
z = y(delay:M:end);  % Downsample with correct offset
plot(1:length(z), z, '.');
xlabel('Symbol Index');
ylabel('Soft Decision Value');
title('Constellation History - 4-PAM Symbols');
grid on;
hold on;
% Draw ideal symbol levels
plot([1 length(z)], [3 3], 'r--');
plot([1 length(z)], [1 1], 'r--');
plot([1 length(z)], [-1 -1], 'r--');
plot([1 length(z)], [-3 -3], 'r--');
legend('Soft Decisions', 'Ideal Levels');
