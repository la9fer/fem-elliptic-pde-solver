% Solves:
%   -Delta u = f       in Omega
%   u = prescribed     on Dirichlet boundary
%   du/dn = g          on Neumann boundary
%
% Uses uniform red refinement and computes L2 and H1-seminorm errors.

clear;
clc;
close all;

format long;

% Generate initial mesh and boundary information.
[Coord, Elem, Nb, Db] = InitialMesh(1);

maxlevel = 5;

L2e = zeros(1, maxlevel);
H1e = zeros(1, maxlevel);
h = zeros(1, maxlevel);

for nl = 1:maxlevel

    % Construct edge-node-element connectivity.
    [n2ed, ed2el] = edge(Elem, Coord);

    % Uniform red refinement.
    [Coord, Elem, Db, Nb] = redrefine( ...
        Coord, Elem, n2ed, ed2el, Db, Nb);

    % Calculate maximum element diameter.
    h(nl) = 0;

    for j = 1:size(Elem, 1)
        vertices = Coord(Elem(j,:), :);

        d12 = norm(vertices(1,:) - vertices(2,:));
        d23 = norm(vertices(2,:) - vertices(3,:));
        d31 = norm(vertices(3,:) - vertices(1,:));

        h(nl) = max(h(nl), max([d12, d23, d31]));
    end

    % Identify free and constrained nodes.
    nNodes = size(Coord, 1);

    FullNodes = 1:nNodes;
    Dbnodes = unique(Db(:));
    FreeNodes = setdiff(FullNodes, Dbnodes);

    % Initialize global stiffness matrix and load vector.
    A = sparse(nNodes, nNodes);
    b = zeros(nNodes, 1);

    % Assemble global stiffness matrix.
    for j = 1:size(Elem, 1)

        nodes = Elem(j, :);

        A(nodes, nodes) = A(nodes, nodes) + ...
            stima(Coord(nodes, :));
    end

    % Assemble volume load vector using centroid quadrature.
    for j = 1:size(Elem, 1)

        nodes = Elem(j, :);
        vertices = Coord(nodes, :);

        area = 0.5 * abs(det([1 1 1; vertices']));
        centroid = sum(vertices, 1) / 3;

        b(nodes) = b(nodes) + ...
            (area / 3) * f(centroid);
    end

    % Assemble Neumann boundary contributions, if present.
    if ~isempty(Nb)

        for j = 1:size(Nb, 1)

            P1 = Coord(Nb(j, 1), :);
            P2 = Coord(Nb(j, 2), :);

            edgeLength = norm(P2 - P1);
            normalDerivative = u_N(P1, P2);

            contribution = edgeLength * normalDerivative / 2;

            b(Nb(j, :)) = b(Nb(j, :)) + contribution;

        end
    end

    % Initialize solution and impose Dirichlet boundary values.
    uh = zeros(nNodes, 1);

    for k = 1:numel(Dbnodes)
        node = Dbnodes(k);
        uh(node) = ue(Coord(node, :));
    end

    % Account for prescribed boundary values.
    b = b - A * uh;

    % Solve the reduced linear system.
    uh(FreeNodes) = A(FreeNodes, FreeNodes) \ b(FreeNodes);

    % Evaluate exact solution at all mesh nodes.
    u = u_nodes(Coord);

    % Compute numerical errors.
    [L2e(nl), H1e(nl)] = Err(Coord, Elem, uh, u);
end

% Calculate observed convergence rates.
ocl2 = zeros(1, maxlevel - 1);
och1 = zeros(1, maxlevel - 1);

for k = 1:maxlevel-1

    ocl2(k) = log(L2e(k) / L2e(k+1)) / ...
              log(h(k) / h(k+1));

    och1(k) = log(H1e(k) / H1e(k+1)) / ...
              log(h(k) / h(k+1));
end

% Display convergence results.
fprintf('\nP1 FEM Convergence Study\n');
fprintf('------------------------------------------------------\n');
fprintf('%-8s %-12s %-14s %-14s\n', ...
    'Level', 'h', 'L2 error', 'H1 error');

for k = 1:maxlevel
    fprintf('%-8d %-12.4e %-14.4e %-14.4e\n', ...
        k, h(k), L2e(k), H1e(k));
end

fprintf('\nObserved convergence rates\n');
fprintf('----------------------------------------\n');

for k = 1:maxlevel-1
    fprintf('Level %d -> %d: L2 = %.4f, H1 = %.4f\n', ...
        k, k+1, ocl2(k), och1(k));
end

fprintf('----------------------------------------\n');

% Display the computed and exact solutions on the final mesh.
show(Coord, Elem, uh, u);
