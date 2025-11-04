clear; clc; 

%-- Define the functions--
% The "@(x)" syntax is a way to create a quick , anonymous fuction 
P= @(x) x.^2 + 2*x -3 ;
P_prime = @(x) 2*x + 2;         % The derivate of P(x)

%------ 2. Find the First Root(by starting with a guess)
x_k =10;
num_iterations = 10;  % repeat the process 10 times 

% fprintf is a formatted print commad
% '%d' is a placeholder for an integer 
fprintf('----Find firt root(starting at x = %d) ---\n', x_k);
for k=1:num_iterations
    %Calculate P(x_k) using our current guess 
    P_val = P(x_k);

    %Calculate P'(x_k) using our current guess 
    P_prime_val = P_prime(x_k);

    %This is the correction term: P(x_k)/ P'(x_k)
    correction = P_val / P_prime_val;

    % This is the full Netwon's formula: x_k+1= x_k-correction 
    x_k_next = x_k - correction ; 

    % Print our prgoress
    % '%2d' means "an integer, 2 spaces wide" (for 'k')
    % '%2d' means "an integer, 2 spaces wide" (for 'k')
    fprintf('Iteration %2d: x = %.6f\n',k,x_k_next);

    %Update the guess to next 
    x_k=x_k_next;
end

