% --- MATLAB/Octave Script for Error Surface: Part (b) ---
% --- Implements J(a) = avg{(y^2[k] - Y^2)^2} ---

clear; clc;

% --- 1. Set up parameters (same as agcerrorsurf.m) ---
n=10000;              % number of steps in simulation
r=randn(n,1);         % generate random inputs (called 'x' in diagram)
ds=0.15;              % desired power of output (Y^2)
range=[-0.7:0.02:0.7]; % range specifies range of values of a
Jagc=zeros(size(range));
j=0;

fprintf('--- Plotting error surface for J(a) = avg{(a^2*r^2 - ds)^2} ---\n');

% --- 2. Run the loop to calculate the surface ---
for a=range           % for each value a
  j=j+1;
  tot=0;
  for i=1:n
    % --- This is the MODIFIED objective function ---
    tot = tot + (a^2 * r(i)^2 - ds)^2;
    % --- End of modification ---
  end
  Jagc(j)=tot/n;      % take average value, and save
end

% --- 3. Plot the results ---
plot(range, Jagc)
ylabel('Cost J(a)')
xlabel('Adaptive gain a')
title('Error Surface for J(a) = avg\{(y^2 - Y^2)^2\}')
grid on;

fprintf('Done.\n');