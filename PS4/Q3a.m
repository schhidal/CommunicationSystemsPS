% --- MATLAB/Octave Script for Steepest Descent: Part (a) ---

clear; clc;

% --- 1. Define the derivative (the slope) ---
% We need the derivative of J(x) = x^2 + 2x - 3
% which is J'(x) = 2x + 2
J_prime = @(x) 2*x + 2;

% --- 2. Set up the iteration parameters ---
N = 100;           % number of iterations
mu = 0.05;         % algorithm step size (a small positive number)
x = zeros(1,N);    % initialize x vector
x(1) = 3;          % starting point x(1) (our initial guess)

fprintf('--- Finding global minimum of J(x) = x^2 + 2x - 3 ---\n');

% --- 3. Run the Steepest Descent loop ---
for k = 1:N-1
    % This is the formula: x_k_next = x_k - mu * J'(x_k)
    x(k+1) = x(k) - mu * J_prime(x(k)); 
end

fprintf('Global Minimum found at: x = %.6f\n', x(N));