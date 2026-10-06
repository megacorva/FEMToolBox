function interpolation=FEInterpolate(mesh, f)
    % returns a FEFunc f_h s.t. f_h=f on the nodes of the mesh
    %-----------------------------------------------------------------
    % inputs:
    %   mesh:P1Mesh
    %   f=@(x,y)...
    %-----------------------------------------------------------------
    % outputs:
    %   interpolation:FEFunc
    %-----------------------------------------------------------------
    nodeValues = arrayfun(f,mesh.nodes(1, :), mesh.nodes(2, :))';
    interpolation = FEFunc(mesh, nodeValues);
end