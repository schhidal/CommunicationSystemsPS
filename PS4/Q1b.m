% --- MATLAB/Octave Script for Newton's Method: Part (b) ---

clear; clc;

% --- 1. Define the functions ---
% We want to find the roots of P'(x), so our functions are:
g       = @(x) 6*x.^2 - 18*x + 12;   % This is P'(x)
g_prime = @(x) 12*x - 18;            % This is P''(x)

% --- 2. Find the First Root (Min/Max) ---
% Let's start with an initial guess
x_k = 10;
num_iterations = 10;

fprintf('--- Finding first min/max (starting at x = %d) ---\n', x_k);
for k = 1:num_iterations
    % Apply Newton's formula
    x_k = x_k - (g(x_k) / g_prime(x_k));
    
    fprintf('Iteration %2d: x = %.6f\n', k, x_k);
end
fprintf('Found Root 1 (a min/max): %.6f\n\n', x_k);


% --- 3. Find the Second Root (Min/Max) ---
% Let's use a different initial guess
x_k = -10;

fprintf('--- Finding second min/max (starting at x = %d) ---\n', x_k);
for k = 1:num_iterations
    % Apply Newton's formula
    x_k = x_k - (g(x_k) / g_prime(x_k));
    
    fprintf('Iteration %2d: x = %.6f\n', k, x_k);
end
fprintf('Found Root 2 (a min/max): %.6f\n', x_k);