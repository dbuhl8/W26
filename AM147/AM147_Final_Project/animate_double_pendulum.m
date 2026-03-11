function animate_double_pendulum(t, y, params, method_name)
% ANIMATE_DOUBLE_PENDULUM  Live animation of double pendulum trajectory
%
%   Inputs:
%     t           : Nx1 time vector
%     y           : Nx4 state matrix [theta1, theta1_dot, theta2, theta2_dot]
%     params      : Struct with fields L1, L2
%     method_name : String label for the integration method (e.g. 'RK4')

L1 = params.L1;
L2 = params.L2;
n  = length(t);

% --- Set up figure ---
figure('Position', [100, 100, 600, 700], 'Color', 'w');

% Axis limits with a little padding
lim = L1 + L2 + 0.1;
set(gca, 'xlim', [-lim lim], 'ylim', [-lim lim*0.4], ...
    'XTick', [-2 -1 0 1 2], 'YTick', [-2 -1 0 1]);
axis square;
grid on;
title(sprintf('Double Pendulum Animation — %s  (L_2 = %.2f m)', ...
      method_name, L2), 'FontSize', 12);
xlabel('x (m)'); ylabel('y (m)');
hold on;

% Pivot point
plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);

% Animated objects
rod1  = animatedline('Color', [0.2157 0.4941 0.7216], ...
                     'LineStyle', '-', 'LineWidth', 3);
bob1  = animatedline('Color', [0.2157 0.4941 0.7216], ...
                     'Marker', 'o', 'MarkerSize', 20, ...
                     'MarkerFaceColor', [0.2157 0.4941 0.7216]);

rod2  = animatedline('Color', [0.8500 0.3250 0.0980], ...
                     'LineStyle', '-', 'LineWidth', 3);
bob2  = animatedline('Color', [0.8500 0.3250 0.0980], ...
                     'Marker', 'o', 'MarkerSize', 20, ...
                     'MarkerFaceColor', [0.8500 0.3250 0.0980]);

% Trailing trace for mass 2 (shows chaotic path)
trace = animatedline('Color', [0.8500 0.3250 0.0980 0.3], ...  % semi-transparent
                     'LineStyle', '-', 'LineWidth', 0.8, ...
                     'MaximumNumPoints', 150);

% Time label
time_text = text(-lim*0.95, lim*0.35, '', 'FontSize', 11, 'FontWeight', 'bold');

legend('', '', '', '', '', ...
       'Trace (m_2)', 'Location', 'northeast');

% --- Animate ---
for k = 1:n
    % Cartesian positions
    x1 =  L1 * sin(y(k,1));
    y1 = -L1 * cos(y(k,1));
    x2 =  x1 + L2 * sin(y(k,3));
    y2 =  y1 - L2 * cos(y(k,3));

    % Update rod 1 (pivot -> mass 1)
    clearpoints(rod1);
    addpoints(rod1, [0, x1], [0, y1]);

    % Update mass 1
    clearpoints(bob1);
    addpoints(bob1, x1, y1);

    % Update rod 2 (mass 1 -> mass 2)
    clearpoints(rod2);
    addpoints(rod2, [x1, x2], [y1, y2]);

    % Update mass 2
    clearpoints(bob2);
    addpoints(bob2, x2, y2);

    % Update trace (keeps last 150 points)
    addpoints(trace, x2, y2);

    % Update time label
    set(time_text, 'String', sprintf('t = %.2f s', t(k)));

    % Render frame
    drawnow;
    pause(t(2) - t(1));   % pause one time step for real-time playback
end

hold off;
end
