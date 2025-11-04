% --- agcvsfading_prob4.m ---
% --- Compares Problem 4's algorithm on the same fading profile ---

clear; clc;

% --- 1. Set up parameters (same as agcvsfading.m) ---
n=50000;              % number of steps in simulation
r=randn(n,1);         % generate raw random inputs
env=0.75+abs(sin(2*pi*[1:n]'/n)); % Fading profile (a)
r=r.*env;             % apply profile to raw input r[k]
ds=0.5;               % desired power of output (Y^2)
a=zeros(1,n); a(1)=1; % initialize AGC parameter
s=zeros(1,n);         % initialize outputs

% --- 2. This is the MODIFIED part ---
mu=0.0001;            % algorithm stepsize (much smaller mu is needed)

fprintf('--- Running AGC for J(a) = avg{(y^2 - Y^2)^2} ---\n');

for k=1:n-1
  s(k)=a(k)*r(k);     % normalize by a(k) to get s[k]
  
  % --- MODIFIED UPDATE RULE from Problem 4 ---
  s_sq = s(k)^2;
  if abs(a(k)) < 1e-4  % Check to avoid divide-by-zero
      update_term = 0;
  else
      % This is the derivative: (y^2/a) * (y^2 - ds)
      update_term = (s_sq / a(k)) * (s_sq - ds);
  end
  a(k+1)=a(k)-mu*update_term;  % Note: No averaging ('avec')
  % --- End of modification ---
end

fprintf('Done.\n');
final_output_power = mean(s(n-1000:n-1).^2);
fprintf('Desired Power (Y^2): %.3f\n', ds);
fprintf('Actual Avg Output Power: %.3f\n', final_output_power);

% --- 3. Plot the results (same as agcvsfading.m) ---
figure;
subplot(3,1,2), plot(a,'g') 
axis([0,length(r),0,1.5])
title('Adaptive gain "a" (Problem 4 Algorithm)')
subplot(3,1,1), plot(r,'r') 
axis([0,length(r),-7,7])
title('Input r(k) - Profile (a)')
subplot(3,1,3),plot(s,'b')
axis([0,length(r),-7,7])
title('Output s(k)')
xlabel('iterations')