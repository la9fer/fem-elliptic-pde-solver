% Convergence study for the P1 finite element method
% applied to the 2D Poisson equation on the unit square.
%
% Computes L2 and H1-seminorm errors over six refinement levels.

clear;
clc;

% Initial mesh and boundary information.
[Coord, Elem, Nb, Db] = InitialMesh(1);

maxlevel = 6;

L2e = zeros(1, maxlevel);
H1e = zeros(1, maxlevel);
h = zeros(1, maxlevel);

for nl = 1:maxlevel

    % Construct edge-node-element connectivity.
    [n2ed, ed2el] = edge(Elem, Coord);

    % Uniform red refinement.
    [Coord, Elem, Db, Nb] = redrefine( ...
        Coord, Elem, n2ed, ed2el, Db, Nb);

    % Estimate mesh size using the area of the first triangle.
    firstElement = Coord(Elem(1,:), :);
    h(nl) = sqrt(abs(det([1 1 1; firstElement']) / 2));

    % Identify free and constrained nodes.
    FullNodes = 1:size(Coord, 1);
    Dbnodes = unique(Db);
    FreeNodes = setdiff(FullNodes, Dbnodes);

    % Initialize global stiffness matrix and load vector.
    nNodes = size(Coord, 1);
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

        area = abs(det([1 1 1; vertices']) / 2);
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
            midpoint = (P1 + P2) / 2;

            normalDerivative = u_N(P1, P2);

            % Midpoint quadrature with linear edge basis functions.
            contribution = edgeLength * normalDerivative / 2;

            b(Nb(j, :)) = b(Nb(j, :)) + contribution;

        end
    end

    % Apply Dirichlet boundary conditions.
    uh = zeros(nNodes, 1);

    if ~isempty(Dbnodes)
        uh(Dbnodes) = ue(Coord(Dbnodes, :)');
    end

    % Modify the load vector to account for prescribed values.
    b = b - A * uh;

    % Solve the reduced linear system.
    uh(FreeNodes) = A(FreeNodes, FreeNodes) \ b(FreeNodes);

    % Evaluate the exact solution at mesh nodes.
    u = u_nodes(Coord);

    % Compute L2 and H1-seminorm errors.
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

% Display convergence rates.
fprintf('\nObserved convergence rates\n');
fprintf('----------------------------------\n');

for k = 1:maxlevel-1
    fprintf('Level %d -> %d: L2 = %.4f, H1 = %.4f\n', ...
        k, k+1, ocl2(k), och1(k));
end

fprintf('----------------------------------\n');
