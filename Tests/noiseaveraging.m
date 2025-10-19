N = 1000; % number of samples
x = randn(1, N); % normally distributed noise,
% mean = 0, var = 1
x = x * sqrt(4); % var = 4
power = sum(x .* x) / N % average power of x
% no abs since x is real