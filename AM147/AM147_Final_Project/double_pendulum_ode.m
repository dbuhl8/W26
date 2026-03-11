function dydt = double_pendulum_ode(t, y, params)
% DOUBLE_PENDULUM_ODE  Right-hand side of double pendulum ODEs
%
%   Inputs:
%     t      : Current time (not used, but required for ODE solvers)
%     y      : 4x1 state vector [theta1; theta1_dot; theta2; theta2_dot]
%     params : Struct with fields g, m1, m2, L1, L2
%
%   Output:
%     dydt   : 4x1 derivative vector [d(theta1)/dt; d^2(theta1)/dt^2; ...]

% Extract parameters
g = params.g;
m1 = params.m1;
m2 = params.m2;
L1 = params.L1;
L2 = params.L2;

% Extract state variables
theta1 = y(1);
theta1_dot = y(2);
theta2 = y(3);
theta2_dot = y(4);

% Compute second derivatives using given equations of motion
% theta1'' equation
num1 = -g * (2*m1 + m2) * sin(theta1) ...
       - m2 * g * sin(theta1 - 2*theta2) ...
       - 2 * sin(theta1 - theta2) * m2 * (theta2_dot^2 * L2 + theta1_dot^2 * L1 * cos(theta1 - theta2));
den1 = L1 * (2*m1 + m2 - m2 * cos(2*theta1 - 2*theta2));
theta1_ddot = num1 / den1;

% theta2'' equation
num2 = 2 * sin(theta1 - theta2) * ...
       (theta1_dot^2 * L1 * (m1 + m2) + ...
        g * (m1 + m2) * cos(theta1) + ...
        theta2_dot^2 * L2 * m2 * cos(theta1 - theta2));
den2 = L2 * (2*m1 + m2 - m2 * cos(2*theta1 - 2*theta2));
theta2_ddot = num2 / den2;

% Construct first-order system: dy/dt = [y2; theta1''; y4; theta2'']
dydt = [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot];

end
