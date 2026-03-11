function [t, y] = trapezoid_method(f, tspan, y0, dt)
% TRAPEZOID_METHOD  Trapezoid (implicit midpoint) method for solving ODEs
%
%   Uses predictor-corrector approach:
%     Predictor: Forward Euler
%     Corrector: Trapezoid rule
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

% Trapezoid time stepping (predictor-corrector)
for i = 1:N-1
    y_curr = y(i, :)';
    t_curr = t(i);
    t_next = t(i+1);
    
    % Predictor: Forward Euler
    f_curr = f(t_curr, y_curr);
    y_pred = y_curr + dt * f_curr;
    
    % Corrector: Trapezoid rule
    f_next = f(t_next, y_pred);
    y(i+1, :) = (y_curr + (dt/2) * (f_curr + f_next))';
end

end
