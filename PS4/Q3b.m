% --- MATLAB/Octave Script for Steepest Descent: Part (b) ---

clear; clc;

% --- 1. Define the derivative (the slope) ---
% We need the derivative of J(x) = 2x^3 - 9x^2 + 12x + 1
% which is J'(x) = 6x^2 - 18x + 12
J_prime = @(x) 6*x.^2 - 18*x + 12;

% --- 2. Set up iteration parameters ---
N = 100;           % number of iterations
mu = 0.01;         % algorithm step size (must be small for this function)

% --- TEST 1: Find the Local Minimum ---
x = zeros(1,N);    % initialize x vector
x(1) = 3;          % starting point x(1) (to the right of the peak at x=1)

fprintf('--- Finding local minimum of J(x) = 2x^3 - 9x^2 + 12x + 1 ---\n');
fprintf('--- Test 1: Starting at x = %d ---\n', x(1));

for k = 1:N-1
    % This is the formula: x_k_next = x_k - mu * J'(x_k)
    x(k+1) = x(k) - mu * J_prime(x(k));
end
fprintf('Found Local Minimum at: x = %.6f\n\n', x(N));


% --- TEST 2: "Finding" the Global Minimum ---
x = zeros(1,N);    % re-initialize x vector
x(1) = 0;          % starting point x(1) (to the left of the peak at x=1)

fprintf('--- Finding global minimum of J(x) = 2x^3 - 9x^2 + 12x + 1 ---\n');
fprintf('--- Test 2: Starting at x = %d ---\n', x(1));

for k = 1:N-1
    x(k+1) = x(k) - mu * J_prime(x(k));
    % Print just a few iterations to show it's "running away"
    if k < 10
        fprintf('Iteration %2d: x = %.6f\n', k, x(k+1));
    end
end
fprintf('...\nFinal value after %d iterations: x = %.6f\n', N, x(N));