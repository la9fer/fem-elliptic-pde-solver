function [Coord, Elem, Db, Nb] = redrefine(Coord, Elem, n2ed, ed2el, Db, Nb)
% REDREFINE Uniform red refinement of a triangular mesh.
%
% Each triangle is divided into four smaller triangles by connecting
% the midpoints of its three edges.
%
% Inputs:
%   Coord : Node coordinates
%   Elem  : Triangle connectivity
%   n2ed  : Node-to-edge connectivity
%   ed2el : Edge-to-element connectivity
%   Db    : Dirichlet boundary edges
%   Nb    : Neumann boundary edges
%
% Outputs:
%   Coord : Refined node coordinates
%   Elem  : Refined triangle connectivity
%   Db    : Refined Dirichlet boundary edges
%   Nb    : Refined Neumann boundary edges

nt = size(Elem, 1);
ne = size(ed2el, 1);

% Create one midpoint node for every mesh edge.
marker = zeros(ne, 1);

for i = 1:ne

    inode = ed2el(i, 1);
    enode = ed2el(i, 2);

    midpoint = (Coord(inode, :) + Coord(enode, :)) / 2;

    Coord(end + 1, :) = midpoint;

    marker(i) = size(Coord, 1);
end

% Refine each triangle into four child triangles.
newElem = zeros(4 * nt, 3);

for i = 1:nt

    ct = Elem(i, :);

    % Retrieve edge indices opposite the corresponding vertices.
    edgeIndices = [
        n2ed(ct(2), ct(3));
        n2ed(ct(3), ct(1));
        n2ed(ct(1), ct(2))
    ];

    m1 = marker(edgeIndices(1));
    m2 = marker(edgeIndices(2));
    m3 = marker(edgeIndices(3));

    % Three corner triangles and one central triangle.
    newElem(4*i - 3, :) = [ct(1), m3, m2];
    newElem(4*i - 2, :) = [ct(2), m1, m3];
    newElem(4*i - 1, :) = [ct(3), m2, m1];
    newElem(4*i, :)     = [m1, m2, m3];

end

Elem = newElem;

% Refine Dirichlet boundary edges.
Db = refineBoundaryEdges(Db, n2ed, marker);

% Refine Neumann boundary edges.
Nb = refineBoundaryEdges(Nb, n2ed, marker);

end

function refinedEdges = refineBoundaryEdges(boundaryEdges, n2ed, marker)
% Split each boundary edge into two edges while preserving orientation.

if isempty(boundaryEdges)
    refinedEdges = boundaryEdges;
    return;
end

nBoundary = size(boundaryEdges, 1);

refinedEdges = zeros(2 * nBoundary, 2);

for i = 1:nBoundary

    node1 = boundaryEdges(i, 1);
    node2 = boundaryEdges(i, 2);

    edgeIndex = n2ed(node1, node2);
    midpointNode = marker(edgeIndex);

    refinedEdges(2*i - 1, :) = [node1, midpointNode];
    refinedEdges(2*i, :)     = [midpointNode, node2];

end
end
