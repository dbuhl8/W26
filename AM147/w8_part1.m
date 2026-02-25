function [c_NE, c_QR, res_NE, res_QR, diff_c] = polyLS_compare(t, y, m)
% polyLS_compare (LEARNER TEMPLATE)
% ------------------------------------------------------------
% Fit a degree-m polynomial to data (t,y) using least squares.
% You must compute the coefficients using:
%   (1) Normal Equations
%   (2) QR Decomposition
% Then compare both approaches and plot the fitted curves.
%
% INPUTS:
%   t : (n x 1) data locations
%   y : (n x 1) observations
%   m : polynomial degree
%
% OUTPUTS:
%   c_NE   : (m+1 x 1) coefficients from normal equations
%   c_QR   : (m+1 x 1) coefficients from QR
%   res_NE : residual norm ||A*c_NE - y||_2
%   res_QR : residual norm ||A*c_QR - y||_2
%   diff_c : norm difference ||c_NE - c_QR||_2
%
% REQUIREMENTS / RESTRICTIONS:
%   - Do NOT use: A\b, pinv, lsqminnorm, polyfit, regress, fit
%   - You MAY use: qr(A,0), matrix multiplications, transpose, norm
%   - You MUST create a plot showing:
%       * data points (scatter)
%       * normal equation fit curve
%       * QR fit curve
% ------------------------------------------------------------

%% Part 1: Build the Vandermonde / design matrix A
% A should have columns: [1, t, t^2, ..., t^m]

n = length(t);                 % number of data points

A = zeros(n, m+1);             % TODO: preallocate A

% TODO: fill A(:,1) with ones
A(:,1) = 1.;

% TODO: fill columns 2..m+1 with powers of t
for i=2:m+1
    A(:,i) = t.^(i-1);
end


%% Part 2: Normal Equations solution
% Compute ATA = A'*A and ATy = A'*y, then solve ATA*c_NE = ATy

% TODO: compute ATA
ATA = A'*A;
% TODO: compute ATy
ATy = A'*y;
% TODO: solve for c_NE (do not use A\b)
c_NE = ATA\ATy;


%% Part 3: QR solution
% Use [Q,R] = qr(A,0), then solve R*c_QR = Q'*y

% TODO: compute [Q,R] = qr(A,0)
[Q, R] = qr(A,0);

% TODO: compute rhs = Q'*y
rhs = Q'*y;

% TODO: solve for c_QR
c_QR = R\rhs;


%% Part 4: Residuals and coefficient difference

% TODO: compute res_NE = norm(A*c_NE - y)
res_NE = sqrt(sum((A*c_NE - y).^2))
% TODO: compute res_QR = norm(A*c_QR - y)
res_QR = sqrt(sum((A*c_QR - y).^2))
% TODO: compute diff_c = norm(c_NE - c_QR)
diff_c = sqrt(sum((c_NE-c_QR).^2))


%% Part 5: REQUIRED Plot
% Plot:
%  - scatter(t,y)
%  - normal equation fit curve evaluated on a dense grid
%  - QR fit curve evaluated on a dense grid

% Dense grid for smooth plotting
tplot = linspace(min(t), max(t), 400)';

% TODO: build Aplot (same structure as A but using tplot)
Aplot = zeros(400, m+1);
Aplot(:,1) = 1;
for i = 2:m+1
    Aplot(:,i) = tplot.^(i-1);
end
% TODO: compute yfit_NE = Aplot*c_NE
yfit_NE = Aplot*c_NE;
% TODO: compute yfit_QR = Aplot*c_QR
yfit_QR = Aplot*c_QR;

% TODO: create the figure and plot:
%   figure
%   scatter(...)
%   hold on
%   plot(tplot, yfit_NE, ...)
%   plot(tplot, yfit_QR, ...)
%   legend(...)
%   grid on, labels, title
figure;
scatter(t, y)
hold on
plot(tplot, yfit_NE,'b-')
plot(tplot, yfit_QR,'r--')
legend('Normal', 'QR')

end
