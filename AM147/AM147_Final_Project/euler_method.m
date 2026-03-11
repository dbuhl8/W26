function [t, y] = euler_method(f, tspan, y0, dt)
% EULER_METHOD  Forward Euler method for solving ODEs
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

% Euler time stepping
for i = 1:N-1
    y_curr = y(i, :)';
    dydt = f(t(i), y_curr);
    y(i+1, :) = (y_curr + dt * dydt)';
end

end
