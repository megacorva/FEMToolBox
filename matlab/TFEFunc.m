classdef TFEFunc<handle
    properties(SetAccess=protected, GetAccess=public)
        mesh P1Mesh
        samplingTime=[]
        nodalValues=[]
        
        % samplingTime:1*M vector,[t_1,...,t_M]
        % nodalValues:N*M vector, u(x_i,t_j)=nodalValues(i,j)
    end

    methods
        function obj=TFEFunc(mesh,samplingTime,nodalValues)
            if size(mesh.nodes,2)~=size(nodalValues,1)
                error('FEFunc:DimensionMismatch', ...
                    'the number of nodes must match the dimension of value vector');
            end
            if size(samplingTime,2)~=size(nodalValues,2)
                error('FEFunc:DimensionMismatch', ...
                    'the number of sampled time must match the dimension of value vector');
            end
            obj.mesh = mesh;
            obj.samplingTime = samplingTime;
            obj.nodalValues = nodalValues;
        end

        function visualize(obj,frameRate,titleFig,path)
            model=createpde();
            geometryFromMesh(model, obj.mesh.nodes, obj.mesh.triangles);
            u=obj.nodalValues;
            uMin = min(u(:));
            uMax = max(u(:));
            xMin = min(obj.mesh.nodes(1,:));
            xMax = max(obj.mesh.nodes(1,:));
            yMin = min(obj.mesh.nodes(2,:));
            yMax = max(obj.mesh.nodes(2,:));
            % 1. Initialize VideoWriter
            if nargin==4
                v = VideoWriter(path, 'MPEG-4');
                v.FrameRate = frameRate;
                open(v);
            end
        
            % 2. Setup Figure window
            fig = figure('Visible','on');
            ax = axes(fig);
            
            for i = 1:length(obj.samplingTime)
            
                pdeplot(ax, model, ...
                    "XYData", u(:,i), ...
                    "ZData", u(:,i));
            
                if nargin >= 3
                    title(ax, titleFig);
                end
            
                text(ax, 0.02, 0.98, ...
                    sprintf('t = %.3f, max = %.5f, min = %.5f', ...
                    obj.samplingTime(i), max(u(:,i)), min(u(:,i))), ...
                    'Units','normalized', ...
                    'VerticalAlignment','top');
            
                view(ax, 2);
            
                xlim(ax,[xMin,xMax]);
                ylim(ax,[yMin,yMax]);
                axis(ax,'equal');
                axis(ax,'manual');
            
                xlabel(ax,'x');
                ylabel(ax,'y');
            
                clim(ax,[uMin,uMax]);
            
                drawnow;
            
                if nargin == 4
                    frame = getframe(fig);
                    writeVideo(v,frame);
                end
            end
        end
    
        function newTFEFunc=uniformSpaceRefine(obj)
            % refine each triangle into 4 congruent triangles, without
            % changing function values
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            n=obj.mesh.nodes;
            tri=obj.mesh.triangles;
            v=obj.nodalValues;
            [newNodes,newValues, newTriangles]=refineTFEFunc4(n,v,tri);
            newMesh=P1Mesh(newNodes,newTriangles);
            newTFEFunc=TFEFunc(newMesh,obj.samplingTime,newValues);
        end

        function newTFEFunc=uniformTimeRefine(obj)
            numNodes=size(obj.mesh.nodes,2);
            numTimeSamples=length(obj.samplingTime);

            % construction of new sampling time and new nodal values
            newSamplingTime=zeros(1,2*numTimeSamples-1);
            newValues=zeros( numNodes,2*numTimeSamples-1 );

            for i=1:numTimeSamples
                newSamplingTime(2*i-1)=obj.samplingTime(i);
                newValues(:,2*i-1)=obj.nodalValues(:,i);
            end

            for i=1:numTimeSamples-1
                newSamplingTime(2*i)=( newSamplingTime(2*i-1)+ ...
                    newSamplingTime(2*i+1) )/2;
                newValues(:,2*i)=( newValues(:,2*i-1)+ ...
                    newValues(:,2*i+1) )/2;
            end
            newTFEFunc=TFEFunc(obj.mesh,newSamplingTime,newValues);
        end
    end
end