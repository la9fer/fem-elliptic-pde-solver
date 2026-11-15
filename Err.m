function [L2e, H1e] = Err(Coord, Elem, uh, u)
% ERR Computes the L2 error and H1-seminorm error for P1 FEM
%
% L2e : Approximation of ||u - uh||_L2(Omega)
% H1e : Approximation of |u - uh|_H1(Omega)
%
% The L2 error uses the existing seven point quadrature rule
% The H1 seminorm error uses a three point degree 2 triangle rule

L2sq = 0;
H1sq = 0;

% Three point quadrature rule on the reference triangle
quadPoints = [
    1/6, 1/6;
    2/3, 1/6;
    1/6, 2/3
];

quadWeights = [1/3; 1/3; 1/3];

for j = 1:size(Elem, 1)

    nodes = Elem(j, :);
    vertices = Coord(nodes, :);

    Cuh = uh(nodes);
    Cu = u(nodes);

    P1 = vertices(1, :);
    P2 = vertices(2, :);
    P3 = vertices(3, :);

    area = 0.5 * abs(det([1 P1; 1 P2; 1 P3]));

    % ---------------------------------------------------------
    % L2 error: retain the original seven-point quadrature
    % ---------------------------------------------------------

    mp12 = 0.5 * (P1 + P2);
    mp13 = 0.5 * (P1 + P3);
    mp23 = 0.5 * (P2 + P3);
    cg = (P1 + P2 + P3) / 3;

    uhmp12 = 0.5 * (Cuh(1) + Cuh(2));
    uhmp13 = 0.5 * (Cuh(1) + Cuh(3));
    uhmp23 = 0.5 * (Cuh(2) + Cuh(3));
    uhcg = sum(Cuh) / 3;

    L2sq = L2sq + area / 60 * ( ...
        8 * (ue(mp12) - uhmp12)^2 + ...
        8 * (ue(mp13) - uhmp13)^2 + ...
        8 * (ue(mp23) - uhmp23)^2 + ...
        3 * (Cu(1) - Cuh(1))^2 + ...
        3 * (Cu(2) - Cuh(2))^2 + ...
        3 * (Cu(3) - Cuh(3))^2 + ...
        27 * (ue(cg) - uhcg)^2 );

    % ---------------------------------------------------------
    % H1 seminorm error: three-point triangle quadrature
    % ---------------------------------------------------------

    A = [ones(1, 3); vertices'];
    basisCoefficients = A' \ eye(3);

    gradients = basisCoefficients(2:3, :);

    % P1 FEM gradient is constant within each triangle
    gradUh = gradients * Cuh;

    for q = 1:3

        xi = quadPoints(q, 1);
        eta = quadPoints(q, 2);

        % Map reference triangle quadrature point to physical triangle
        point = (1 - xi - eta) * P1 + ...
                xi * P2 + eta * P3;

        [uxExact, uyExact] = uxe(point);

        gradExact = [uxExact; uyExact];

        gradError = gradExact - gradUh;

        H1sq = H1sq + area * quadWeights(q) * ...
               sum(gradError.^2);

    end
end

L2e = sqrt(max(L2sq, 0));
H1e = sqrt(max(H1sq, 0));

end
