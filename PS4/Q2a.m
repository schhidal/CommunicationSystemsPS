% --- MATLAB/Octave Script for Steepest Descent: Part (a) ---

clear; clc;

% --- 1. Define the derivative (the slope) ---
% We need the derivative of P(x) = x^2 + 2x - 3
% which is P'(x) = 2x + 2
P_prime = @(x) 2*x + 2;

% --- 2. Set up the iteration parameters ---
N = 100;           % number of iterations
mu = 0.05;         % algorithm step size (you can try 0.01 or 0.1 too)
x = zeros(1,N);    % initialize x to zero (like in polyconverge.m)
x(1) = 10;         % starting point x(1) (our initial guess)

fprintf('--- Finding minimum of P(x) = x^2 + 2x - 3 ---\n');
fprintf('--- Starting at x = %d, step size mu = %.2f ---\n', x(1), mu);

% --- 3. Run the Steepest Descent loop ---
for k = 1:N-1
    % This is the general formula from the problem:
    % x_k_next = x_k - mu * P'(x_k)
    x(k+1) = x(k) - mu * P_prime(x(k));
    
    % This line is just to stop printing if it converges early
    if abs(x(k+1) - x(k)) < 1e-6
        fprintf('Converged at iteration %d\n', k);
        x(k+2:end) = x(k+1); % Fill rest of array with final value
        break; % Stop the loop
    end
end

fprintf('Final Minimum found at: x = %.6f\n', x(N));