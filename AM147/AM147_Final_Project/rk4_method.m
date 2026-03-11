function [t, y] = rk4_method(f, tspan, y0, dt)
% RK4_METHOD  Classical 4th-order Runge-Kutta method for solving ODEs
%
%   Inputs:
%     f      : Function handle for dy/dt = f(t, y)
%     tspan  : [t0, tf] time interval
%     y0     : Initial condition (column vector)
%     dt     : Time step size
%
%   Outputs:
%     t : Column vector of time points
%     y : Matrix where row i is the solution at t(i)

% Initialize time vector
t = (tspan(1):dt:tspan(2))';
N = length(t);

% Initialize solution matrix
n = length(y0);
y = zeros(N, n);
y(1, :) = y0';

% RK4 time stepping
for i = 1:N-1
    y_curr = y(i, :)';
    t_curr = t(i);
    
    % Compute the four stages
    s1 = f(t_curr, y_curr);
    s2 = f(t_curr + dt/2, y_curr + (dt/2)*s1);
    s3 = f(t_curr + dt/2, y_curr + (dt/2)*s2);
    s4 = f(t_curr + dt, y_curr + dt*s3);
    
    % Weighted average
    y(i+1, :) = (y_curr + (dt/6) * (s1 + 2*s2 + 2*s3 + s4))';
end

end
