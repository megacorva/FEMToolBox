function [newNodes,newValues,newTriangles] = refineTFEFunc4(nodes,values,triangles)
% Uniformly refine a triangular mesh, keeping the FE function value:
% each triangle is divided into 4 congruent/similar triangles.
%
% Input:
%   nodes       : 2 x Nn
%   values      : Nn x 1
%   triangles   : 3 x Nt
%
% Output:
%   newNodes     : 2 x newNn
%   newTriangles : 3 x (4*Nt)

    Nn = size(nodes,2);
    Nt = size(triangles,2);

    % -----------------------------------------
    % 1. Collect all edges of all triangles
    % -----------------------------------------
    edges = [triangles([1,2],:), ...
             triangles([2,3],:), ...
             triangles([3,1],:)];

    % Represent each edge canonically as (min,max)
    edges = sort(edges,1);

    % -----------------------------------------
    % 2. Find unique edges
    % -----------------------------------------
    [uniqueEdges,~,edgeID] = unique(edges','rows');
    uniqueEdges = uniqueEdges';

    Ne = size(uniqueEdges,2);

    % -----------------------------------------
    % 3. Create midpoint nodes
    % -----------------------------------------
    midpointNodes = ...
        (nodes(:,uniqueEdges(1,:)) + ...
         nodes(:,uniqueEdges(2,:))) / 2;

    midpointValues=(values(uniqueEdges(1,:),:) + ...
         values(uniqueEdges(2,:),:)) / 2;

    newNodes = [nodes, midpointNodes];
    newValues= [values; midpointValues];
    % midpoint node and value corresponding to unique edge k
    midpointIndex = Nn + (1:Ne);

    % edgeID contains:
    % first Nt entries   -> edge (1,2)
    % next Nt entries    -> edge (2,3)
    % last Nt entries    -> edge (3,1)

    m12 = midpointIndex(edgeID(1:Nt));
    m23 = midpointIndex(edgeID(Nt+1:2*Nt));
    m31 = midpointIndex(edgeID(2*Nt+1:3*Nt));

    % -----------------------------------------
    % 4. Split every triangle into four
    % -----------------------------------------
    n1 = triangles(1,:);
    n2 = triangles(2,:);
    n3 = triangles(3,:);

    T1 = [n1;   m12; m31];
    T2 = [m12; n2;   m23];
    T3 = [m31; m23; n3];
    T4 = [m12; m23; m31];

    newTriangles = [T1,T2,T3,T4];
end