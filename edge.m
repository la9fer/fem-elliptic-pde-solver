function [nodes2edge, edge2element] = edge(element, coordinate)
% EDGE Construct node-edge and edge-element connectivity.
%
% Inputs:
%   element    : M-by-3 triangular element connectivity
%   coordinate : N-by-2 node coordinates
%
% Outputs:
%   nodes2edge : Sparse symmetric node-to-edge connectivity matrix
%   edge2element : E-by-4 matrix containing:
%       columns 1-2: edge endpoint nodes
%       columns 3-4: adjacent element information

nNodes = size(coordinate, 1);
nElements = size(element, 1);

% Collect all triangle edges.
allEdges = [
    element(:, [1 2]);
    element(:, [2 3]);
    element(:, [3 1])
];

% Sort endpoints so that each undirected edge has one representation.
sortedEdges = sort(allEdges, 2);

% Identify unique edges.
[uniqueEdges, ~, edgeIndex] = unique(sortedEdges, 'rows');

nEdges = size(uniqueEdges, 1);

% Construct symmetric node-to-edge connectivity.
nodes2edge = sparse( ...
    [uniqueEdges(:,1); uniqueEdges(:,2)], ...
    [uniqueEdges(:,2); uniqueEdges(:,1)], ...
    [1:nEdges, 1:nEdges], ...
    nNodes, nNodes);

% Construct edge-to-element connectivity.
edge2element = zeros(nEdges, 4);

% Element indices associated with each edge.
elementIndices = repmat((1:nElements)', 3, 1);

for k = 1:size(allEdges, 1)

    e = edgeIndex(k);
    elemIndex = elementIndices(k);

    if edge2element(e, 1) == 0

        edge2element(e, 1:2) = allEdges(k, :);
        edge2element(e, 3) = elemIndex;

    else

        edge2element(e, 4) = elemIndex;

    end
end
end
