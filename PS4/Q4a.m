% --- MATLAB/Octave Script for Steepest Descent: Part (a) ---
% --- Implements J(a) = avg{(y^2[k] - Y^2)^2} ---

clear; clc;

% --- 1. Set up parameters (same as agcgrad.m) ---
n = 10000;              % number of steps in simulation
vr = 1.0;               % power (variance) of the input
x = sqrt(vr)*randn(n,1); % input signal (I use 'x' to match diagram)
Y_sq = 0.15;            % desired output power (Y^2), called 'ds' in agcgrad
mu = 0.0001;            % algorithm stepsize (needs to be smaller for this alg)
lenavg = 10;            % length over which to average
a = zeros(n,1); a(1)=1; % initialize AGC parameter
y = zeros(n,1);         % initialize outputs
avec = zeros(1,lenavg); % vector to store terms for averaging

fprintf('--- Running AGC for J(a) = avg{(y^2 - Y^2)^2} ---\n');

% --- 2. Run the Steepest Descent loop ---
for k = 1:n-1
    % Calculate the output
    y(k) = a(k) * x(k);
    
    % --- This is our NEW update term ---
    y_sq = y(k)^2;
    
    % Check for a(k) being too small, to avoid divide-by-zero
    if abs(a(k)) < 1e-4 
        update_term = 0;
    else
        % This is the derivative we derived: (y^2/a) * (y^2 - Y^2)
        update_term = (y_sq / a(k)) * (y_sq - Y_sq);
    end
    % --- End of new part ---

    % This is the averaging from your agcgrad.m
    avec = [update_term, avec(1:lenavg-1)]; 
    
    % This is the update from your agcgrad.m
    a(k+1) = a(k) - mu * mean(avec); 
end

fprintf('Done.\n');
% Let's also check the final average power of the output
final_output_power = mean(y(n-1000:n-1).^2);
fprintf('Desired Power (Y^2): %.3f\n', Y_sq);
fprintf('Actual Avg Output Power: %.3f\n', final_output_power);

% --- 3. Plot the results (same as agcgrad.m) ---
figure;
subplot(3,1,1)
plot(a)
title('Adaptive gain "a" (for J = avg{(y^2-Y^2)^2})')
subplot(3,1,2)
plot(x,'r') 
axis([0,n,-5,5])
title('Input x(k)')
subplot(3,1,3)
plot(y,'b')
axis([0,n,-5,5])
title('Output y(k)')
xlabel('iterations')