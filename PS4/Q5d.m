% --- agcvsfading_default_d.m ---
% --- Tests the "default" algorithm on fading profile (d) ---

clear; clc;

% Define a simple ustep function
ustep = @(x) (x >= 0);

% --- 1. Set up parameters (same as agcvsfading.m) ---
n=50000;
r=randn(n,1);

% --- THIS IS THE MODIFIED LINE FOR PART (d) ---
k_vec = [1:n]';
env = ustep(k_vec) - ustep(k_vec-10000) + ustep(k_vec-20000);
% --- End of modification ---

r=r.*env;
ds=0.5;
a=zeros(1,n); a(1)=1;
s=zeros(1,n);
mu=0.01;

% --- 2. Run the "Default" loop ---
for k=1:n-1
  s(k)=a(k)*r(k);
  a(k+1)=a(k)-mu*(s(k)^2-ds); % Default update rule
end

% --- 3. Plot the results (same as agcvsfading.m) ---
figure;
subplot(3,1,2), plot(a,'g') 
title('Adaptive gain "a" (Default Algorithm)')
subplot(3,1,1), plot(r,'r') 
axis([0,length(r),-7,7])
title('Input r(k) - Profile (d)')
subplot(3,1,3),plot(s,'b')
axis([0,length(r),-7,7])
title('Output s(k)')
xlabel('iterations')