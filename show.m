function show(coordinates, elements, uh, u)
% SHOW Visualize the finite-element and exact solutions
%
% Inputs:
%   coordinates : N-by-2 matrix of node coordinates
%   elements    : M-by-3 matrix of triangular element connectivity
%   uh          : finite-element solution at the nodes
%   u           : exact solution at the nodes

    uh = uh(:);
    u  = u(:);

    % Use a common color scale for a fair visual comparison
    allValues = [uh; u];
    colorLimits = [min(allValues), max(allValues)];

    % Avoid identical color limits for a constant solution.
    if colorLimits(1) == colorLimits(2)
        colorLimits = colorLimits + [-1, 1] * eps(max(1, abs(colorLimits(1))));
    end

    figure('Color', 'w', 'Name', 'FEM Solution Comparison');

    % Finite-element solution.
    subplot(1, 2, 1);
    trisurf(elements, coordinates(:,1), coordinates(:,2), uh, ...
        'FaceColor', 'interp', ...
        'EdgeColor', [0.35 0.35 0.35]);

    view(2);
    axis equal tight;
    caxis(colorLimits);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Finite Element Solution');
    set(gca, 'YDir', 'normal');

    % Exact solution.
    subplot(1, 2, 2);
    trisurf(elements, coordinates(:,1), coordinates(:,2), u, ...
        'FaceColor', 'interp', ...
        'EdgeColor', [0.35 0.35 0.35]);

    view(2);
    axis equal tight;
    caxis(colorLimits);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Exact Solution');
    set(gca, 'YDir', 'normal');

    colormap(parula);
end
