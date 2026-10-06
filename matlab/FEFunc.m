classdef FEFunc<handle
    properties(SetAccess=protected, GetAccess=public)
        mesh P1Mesh
        nodalValues=[]

        % nodalValues:N*1 vector, i-th component is taken at the i-th
        % node
    end

    methods
        function obj=FEFunc(mesh,nodalValues)
            if size(mesh.nodes,2)~=size(nodalValues,1)
                error('FEFunc:DimensionMismatch', ...
                    'the number of nodes must match the dimension of value vector');
            end
            obj.mesh = mesh;
            obj.nodalValues = nodalValues;
        end

        function visualize(obj,titleFig,path)
            model=createpde();
            geometryFromMesh(model, obj.mesh.nodes, obj.mesh.triangles);
            pdeplot(model, 'XYData', obj.nodalValues);
            if nargin>=2
                title(titleFig);
            end
            xlabel("x");
            ylabel("y");
            axis equal;
            
            % Save only if path is provided
            if nargin == 3
                exportgraphics(gcf, path, "Resolution", 300);
            end
        end
    
        function newFEFunc=uniformRefine(obj)
            % refine each triangle into 4 congruent triangles, without
            % changing function values
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            n=obj.mesh.nodes;
            t=obj.mesh.triangles;
            v=obj.nodalValues;
            [newNodes,newValues, newTriangles]=refineFEFunc4(n,v,t);
            newMesh=P1Mesh(newNodes,newTriangles);
            newFEFunc=FEFunc(newMesh,newValues);
        end
    end
end