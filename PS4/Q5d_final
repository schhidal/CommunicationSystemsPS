% --- agcvsfading_prob4_d.m ---
% --- Tests the "Problem 4" algorithm on fading profile (d) ---

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
mu=0.0001; % Using the smaller mu for this algorithm

% --- 2. Run the "Problem 4" loop ---
for k=1:n-1
  s(k)=a(k)*r(k);
  
  % --- MODIFIED UPDATE RULE from Problem 4 ---
  s_sq = s(k)^2;
  if abs(a(k)) < 1e-4
      update_term = 0;
  else
      update_term = (s_sq / a(k)) * (s_sq - ds);
  end
  a(k+1)=a(k)-mu*update_term;
  % --- End of modification ---
end

% --- 3. Plot the results (same as agcvsfading.m) ---
figure;
subplot(3,1,2), plot(a,'g') 
title('Adaptive gain "a" (Problem 4 Algorithm)')
subplot(3,1,1), plot(r,'r') 
axis([0,length(r),-7,7])
title('Input r(k) - Profile (d)')
subplot(3,1,3),plot(s,'b')
axis([0,length(r),-7,7])
title('Output s(k)')
xlabel('iterations')