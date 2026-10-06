function newMesh=editP1Mesh(mesh)
    % edits the given mesh, and returns the edited new P1Mesh object
    % inputs and outputs---------------------------------------------------
    %   if the input is empty, then the editor starts from a blank mesh,
    %   both a .mat file and a P1Mesh can be used to call the python
    %   editor.
    projectRoot = fullfile(fileparts(mfilename('fullpath')), '..');
    if count(py.sys.path, projectRoot) == 0
        insert(py.sys.path, int32(0), projectRoot);
    end
    editor = py.importlib.import_module('editP1Mesh');
    
    if nargin == 0
        nodes = [];
        triangles = [];

    elseif isstring(mesh) || ischar(mesh)
        data = load(mesh);
        nodes = data.nodes;
        triangles = data.triangles;
    
    elseif isa(mesh, 'P1Mesh')
        nodes = mesh.nodes;
        triangles = mesh.triangles;
    
    else
        error('editP1Mesh:InvalidInput', ...
            'Input must be a MAT-file path or a P1Mesh object.');
    end

    res=editor.editP1Mesh( py.numpy.array(nodes), ...
            py.numpy.array(triangles) );
    newNodes=double(res{1});
    newTriangles=double(res{2});
    newMesh=P1Mesh(newNodes,newTriangles);
end