%% AM 147 Final Project - Double Pendulum Simulation
% Simulates a double pendulum using Euler, Trapezoid, and RK4 methods
% and compares each against MATLAB's ode45 solver on a uniform output grid.

%% Clear Workspace
clear; close all; clc;

%% Parameters

g  = 9.81;                    % gravity (m/s^2)
m1 = 1;  m2 = 1;              % masses (kg)
L1 = 1;                       % rod 1 length (m)
L2_values = [0.25, 0.5];      % rod 2 lengths to test (m)
y0 = [pi/3; 0; 2*pi/3; 0];    % [theta1; theta1_dot; theta2; theta2_dot]
tspan     = [0, 5];           % time interval (s)
dt_values = [0.1, 0.01];      % time steps to test (s)

method_funcs = {@euler_method, @trapezoid_method, @rk4_method};
method_names = {'Euler', 'Trapezoid', 'RK4'};

%% TASK 1: ANIMATE THE DOUBLE PENDULUM

% Use smallest time step and largest L2 for the animation
dt_anim = min(dt_values);
L2_anim = max(L2_values);

% Store system parameters in a structure
params.g  = g;
params.m1 = m1;
params.m2 = m2;
params.L1 = L1;
params.L2 = L2_anim;

tgrid_anim = tspan(1):dt_anim:tspan(2);

% ode45 returns the solution evaluated on the uniform output grid
[t_ref, y_ref] = ode45(@(t,y) double_pendulum_ode(t, y, params), tgrid_anim, y0);

animate_double_pendulum(t_ref, y_ref, params, 'ode45 (reference)');

%% TASKS 2-4: METHOD COMPARISON vs ode45 %%

dt_plot    = 0.01;
tgrid_plot = tspan(1):dt_plot:tspan(2);

state_names = {'$\theta_1$', '$\dot{\theta}_1$', '$\theta_2$', '$\dot{\theta}_2$'};
state_units = {'(rad)', '(rad/s)', '(rad)', '(rad/s)'};

clr = [0.8500 0.3250 0.0980;   % Euler     - orange
       0.9290 0.6940 0.1250;   % Trapezoid - yellow
       0.4940 0.1840 0.5560];  % RK4       - purple

for i = 1:length(L2_values)
    L2 = L2_values(i);
    params.L2 = L2;

    % Reference solution sampled on the same uniform grid
    [t_ode45, y_ode45] = ode45(@(t,y) double_pendulum_ode(t, y, params), tgrid_plot, y0);

    % Fixed-step methods
    [t_euler, y_euler] = euler_method(@(t,y) double_pendulum_ode(t, y, params), tspan, y0, dt_plot);
    [t_trap,  y_trap]  = trapezoid_method(@(t,y) double_pendulum_ode(t, y, params), tspan, y0, dt_plot);
    [t_rk4,   y_rk4]   = rk4_method(@(t,y) double_pendulum_ode(t, y, params), tspan, y0, dt_plot);

    % Direct error computation on the common grid
    err_euler = abs(y_euler - y_ode45);
    err_trap  = abs(y_trap  - y_ode45);
    err_rk4   = abs(y_rk4   - y_ode45);

    figure('Position', [100, 100, 1400, 600]);

    for s = 1:4
        subplot(2, 4, s);
        plot(t_euler, y_euler(:,s), '-', 'Color', clr(1,:), 'LineWidth', 1.5); hold on;
        plot(t_trap,  y_trap(:,s),  '-', 'Color', clr(2,:), 'LineWidth', 1.5);
        plot(t_rk4,   y_rk4(:,s),   '-', 'Color', clr(3,:), 'LineWidth', 1.5);
        plot(t_ode45, y_ode45(:,s), 'k--', 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel(sprintf('%s %s', state_names{s}, state_units{s}), 'Interpreter', 'latex');
        title(sprintf('State: %s', state_names{s}), 'Interpreter', 'latex');
        if s == 1
            legend('Euler', 'Trapezoid', 'RK4', 'ode45', 'Location', 'best');
        end
        grid on;

        subplot(2, 4, 4 + s);
        semilogy(t_ode45, err_euler(:,s), '-', 'Color', clr(1,:), 'LineWidth', 1.5); hold on;
        semilogy(t_ode45, err_trap(:,s),  '-', 'Color', clr(2,:), 'LineWidth', 1.5);
        semilogy(t_ode45, err_rk4(:,s),   '-', 'Color', clr(3,:), 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel(sprintf('|Error| in %s', state_names{s}), 'Interpreter', 'latex');
        title(sprintf('Error: %s', state_names{s}), 'Interpreter', 'latex');
        if s == 1
            legend('Euler', 'Trapezoid', 'RK4', 'Location', 'best');
        end
        grid on;
    end

    sgtitle(sprintf('Method Comparison: L_2 = %.2f m, \\Delta t = %.2f s', L2, dt_plot), ...
            'FontSize', 14, 'FontWeight', 'bold');

    saveas(gcf, sprintf('method_comparison_L2_%d.png', round(L2*100)));
end

%% TASK 5: RMSE AND COMPUTATIONAL COST %%

rmse_all  = zeros(3, length(L2_values), length(dt_values), 4);
times_all = zeros(3, length(L2_values), length(dt_values));

for i = 1:length(L2_values)
    L2 = L2_values(i);
    params.L2 = L2;

    for j = 1:length(dt_values)
        dt = dt_values(j);
        tgrid = tspan(1):dt:tspan(2);

        % Reference solution sampled on the same grid
        [t_ref, y_ref] = ode45(@(t,y) double_pendulum_ode(t, y, params), tgrid, y0);

        for k = 1:3
            tic;
            [t_m, y_m] = method_funcs{k}( ...
                @(t,y) double_pendulum_ode(t, y, params), ...
                tspan, y0, dt);
            times_all(k, i, j) = toc;

            % Compute RMSE for each state variable
            for s = 1:4
                rmse_all(k, i, j, s) = sqrt(mean((y_m(:,s) - y_ref(:,s)).^2));
            end
        end
    end
end

%% SUMMARY TABLE %%
fprintf('\n========================================\n');
fprintf('SUMMARY: RMSE (dt = 0.01 s, L2 = 0.50 m)\n');
fprintf('========================================\n\n');

fprintf('%-12s | %-12s %-12s %-12s %-12s | Time (s)\n', ...
        'Method', 'theta1', 'theta1_dot', 'theta2', 'theta2_dot');
fprintf('%s\n', repmat('-', 1, 80));

for k = 1:3
    r = squeeze(rmse_all(k, 2, 2, :))';
    fprintf('%-12s | %-12.4e %-12.4e %-12.4e %-12.4e | %.4f\n', ...
            method_names{k}, r(1), r(2), r(3), r(4), times_all(k, 2, 2));
end
