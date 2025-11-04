% --- MATLAB/Octave Script for Steepest Descent: Part (c) ---

clear; clc;

% --- 1. Define the derivative (the slope) ---
% We need the derivative of J(x) = 1 + 2|x| + 6sin(6x)
% which is J'(x) = 2*sign(x) + 36*cos(6x)
% Note: sign(x) is a built-in MATLAB/Octave function.
J_prime = @(x) 2*sign(x) + 36*cos(6*x);

% --- 2. Set up iteration parameters ---
N = 100;           % number of iterations
mu = 0.001;        % algorithm step size (must be very small!)

% --- TEST 1: Start with a guess near 0 ---
x = zeros(1,N);    % initialize x vector
x(1) = 0.1;        % starting point x(1)
fprintf('--- Finding a local minimum of J(x) = 1 + 2|x| + 6sin(6x) ---\n');
fprintf('--- Test 1: Starting at x = %.1f ---\n', x(1));

for k = 1:N-1
    x(k+1) = x(k) - mu * J_prime(x(k)); 
end
fprintf('Found a Local Minimum at: x = %.6f\n\n', x(N));
% Let's also calculate the value of J(x) at this minimum
J_val_1 = 1 + 2*abs(x(N)) + 6*sin(6*x(N));
fprintf('Value of J(x) at this minimum: J(%.6f) = %.6f\n\n', x(N), J_val_1);


% --- TEST 2: Start with a guess further away ---
x = zeros(1,N);    % re-initialize x vector
x(1) = 0.7;        % new starting point x(1)

fprintf('--- Finding a local minimum of J(x) = 1 + 2|x| + 6sin(6x) ---\n');
fprintf('--- Test 2: Starting at x = %.1f ---\n', x(1));

for k = 1:N-1
    x(k+1) = x(k) - mu * J_prime(x(k));
end
fprintf('Found a Local Minimum at: x = %.6f\n\n', x(N));
% Let's calculate the value of J(x) at this minimum
J_val_2 = 1 + 2*abs(x(N)) + 6*sin(6*x(N));
fprintf('Value of J(x) at this minimum: J(%.6f) = %.6f\n\n', x(N), J_val_2);