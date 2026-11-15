% Runs the P1 FEM convergence experiment, prints the error
% and observed convergence rates, and saves a log-log plot.
%
% Run from the repository directory:
%   convergence_study

clear;
clc;
close all;

% Run the FEM convergence driver.
% This defines L2e, H1e, ocl2 and och1.
ocfem;

% Number of refinement levels.
nlevels = numel(L2e);

% Nominal mesh-size sequence.
h = 1 ./ (2 .^ (1:nlevels));

% Validate the convergence data.
if numel(H1e) ~= nlevels
    error('L2e and H1e must contain the same number of entries.');
end

if numel(ocl2) ~= nlevels - 1 || numel(och1) ~= nlevels - 1
    error('Convergence-rate arrays have unexpected dimensions.');
end

% Print convergence table.
fprintf('\nP1 FEM Convergence Study\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%-7s %-10s %-14s %-14s %-12s %-12s\n', ...
    'Level', '1/h', 'L2 error', 'H1 error', 'L2 rate', 'H1 rate');
fprintf('---------------------------------------------------------------\n');

for k = 1:nlevels

    if k == 1
        fprintf('%-7d %-10d %-14.4e %-14.4e %-12s %-12s\n', ...
            k, round(1/h(k)), L2e(k), H1e(k), '-', '-');
    else
        fprintf('%-7d %-10d %-14.4e %-14.4e %-12.4f %-12.4f\n', ...
            k, round(1/h(k)), L2e(k), H1e(k), ...
            ocl2(k-1), och1(k-1));
    end
end

fprintf('---------------------------------------------------------------\n');

% Plot numerical errors and theoretical reference slopes.
figure;

loglog(h, L2e, 'o-', 'LineWidth', 1.5, ...
    'MarkerSize', 7);
hold on;

loglog(h, H1e, 's-', 'LineWidth', 1.5, ...
    'MarkerSize', 7);

% Reference slopes, anchored at the first error values.
loglog(h, L2e(1) * (h/h(1)).^2, 'k--', ...
    'LineWidth', 1.2);

loglog(h, H1e(1) * (h/h(1)), 'k:', ...
    'LineWidth', 1.2);

xlabel('Nominal mesh size h');
ylabel('Error');
title('P1 FEM Convergence for the 2D Poisson Equation');

legend('L^2 error', 'H^1 seminorm error', ...
    'Reference slope 2', 'Reference slope 1', ...
    'Location', 'southeast');

grid on;
box on;

% Save convergence plot.
print('-dpng', 'convergence.png');

fprintf('\nConvergence plot saved as convergence.png\n');
