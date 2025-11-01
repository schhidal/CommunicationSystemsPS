% --- Part (b): Simulate Step Response using conv ---

clear; clc; close all;

% Define the impulse response h[n] = delta[n-1] - delta[n-2]
% Index: n=0  n=1  n=2
h = [0,   1,  -1];

% Define the unit step input x[n] = u[n] for n=0 to 9
N_x = 10; % Length of input signal
n_x = 0:(N_x-1);
x = ones(1, N_x); % u[n] = 1 for n>=0

% Compute the convolution using the conv function
y = conv(x, h);

% Determine the time index vector for the output y[n]
% The output starts at n = (start index of x) + (start index of h) = 0 + 0 = 0
N_y = length(y); % Length of output is length(x)+length(h)-1
n_y = 0:(N_y-1);

% Plot the input, impulse response, and output
figure('Name', 'Step Response Simulation');

subplot(3, 1, 1);
stem(n_x, x, 'filled');
title('Input x[n] = u[n]');
xlabel('n');
ylabel('Amplitude');
xlim([-1, N_y]); % Adjust x-axis for comparison
grid on;

subplot(3, 1, 2);
stem(0:(length(h)-1), h, 'filled'); % Plot h starting at n=0
title('Impulse Response h[n] = \delta[n-1] - \delta[n-2]');
xlabel('n');
ylabel('Amplitude');
xlim([-1, N_y]);
grid on;

subplot(3, 1, 3);
stem(n_y, y, 'filled');
title('Output y[n] = (x * h)[n]');
xlabel('n');
ylabel('Amplitude');
xlim([-1, N_y]);
grid on;

% Display the numerical output vector y
disp('Output vector y[n]:');
disp(y);