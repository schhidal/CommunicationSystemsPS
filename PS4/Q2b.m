% --- MATLAB/Octave Script for Steepest Descent: Part (b) ---

clear; clc;

% --- 1. Define the derivative (the slope) ---
% We need the derivative of P(x) = 2x^3 - 9x^2 + 12x + 1
% which is P'(x) = 6x^2 - 18x + 12
P_prime = @(x) 6*x.^2 - 18*x + 12;

% --- 2. Set up the iteration parameters ---
N = 100;           % number of iterations
mu = 0.01;         % algorithm step size (this function is steeper, so we use a smaller mu)
x = zeros(1,N);    % initialize x to zero
x(1) = 10;         % starting point x(1) (our initial guess)

fprintf('--- Finding local minimum of P(x) = 2x^3 - 9x^2 + 12x + 1 ---\n');
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