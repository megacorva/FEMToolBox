function boundaryNodes=findBoundaryNodes(nodes,triangles)
    % inputs:
    %    nodes: 2*N matrix
        %[X1,X2,...,XN;
        %[Y1,Y2,...,YN]

        %   triangles: 3*N matrix representing nodes of triangles
        %[n11, n21,......,nN1;
        % n12,...............;
        % n13,...........,nN3]

    % outputs:
    % a 1*n matrix of boundary nodes, sorted in 
    % counterclockwise direction

    % Collect all triangle edges
    edges = [
        triangles([1 2], :)';
        triangles([2 3], :)';
        triangles([3 1], :)'
    ];

    edges = sort(edges,2);
    
    % an edge is on the boundary iff it is adjacent to only 1 triangle
    [uniqueEdges,~,ic] = unique(edges,'rows');

    counts = accumarray(ic,1);
    
    boundaryEdges = uniqueEdges(counts == 1,:);
    

    % Construct adjacency information
    N = size(nodes,2);
    adj = cell(N,1);

    for k = 1:size(boundaryEdges,1)
        i = boundaryEdges(k,1);
        j = boundaryEdges(k,2);

        adj{i}(end+1) = j;
        adj{j}(end+1) = i;
    end

    % Walk around boundary
    numBNodes=size(boundaryEdges,1);
    boundaryNodes = zeros(numBNodes,1);


    % Start with an arbitrary boundary edge
    n1 = boundaryEdges(1,1);
    n2 = boundaryEdges(1,2);

    boundaryNodes(1) = n1;
    boundaryNodes(2) = n2;

    for k = 3:numBNodes
        previous = boundaryNodes(k-2);
        current  = boundaryNodes(k-1);

        neighbors = adj{current};

        % Choose the neighbor that we did not just come from
        next = neighbors(neighbors ~= previous);

        boundaryNodes(k) = next;
    end

    % Check orientation using signed polygon area
    x = nodes(1,boundaryNodes);
    y = nodes(2,boundaryNodes);

    signedArea = sum( ...
        x .* y([2:end 1]) ...
        - y .* x([2:end 1]) );

    % Negative means clockwise
    if signedArea < 0
        boundaryNodes = flip(boundaryNodes);
    end
end