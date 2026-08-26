classdef P1Mesh<handle
    
    properties (SetAccess = protected, GetAccess = public)
        % independent properties:
        nodes=[]
        triangles=[]

        %   nodes: 2*N matrix
        %[X1,X2,...,XN;
        %[Y1,Y2,...,YN]

        %   triangles: 3*N matrix representing nodes of triangles
        %[n11, n21,......,nN1;
        % n12,...............;
        % n13,...........,nN3]
        
        % dependent properties:
        boundaryNodes=[]
        % boundaryNodes: 1*N matrix of node indexes, arranged
        % in the positive (counterclockwise) direction

        internalNodes=[]
        % internalNodes: 1*N matrix of node indexes

        isConvex=true
    end

    methods

        function obj=P1Mesh(nodes,triangles)
            obj.nodes=nodes;
            obj.triangles=triangles;

            % find boundary nodes
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            bNodes=findBoundaryNodes(nodes,triangles);
            obj.boundaryNodes=bNodes;
            numNodes = size(nodes, 2);
            obj.internalNodes=setdiff(1:numNodes,bNodes);

            boundary=nodes(:,bNodes);
            obj.isConvex=testConvexity(boundary);
        end

        function e21=getE21Cons(obj)
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            %|u-Iu|_1<=E21 |u|_2, computed with Kobayashi formulae
            n=obj.nodes;
            t=obj.triangles;
            e21=0;
            for i=1:size(t,2)
                nodeIndeces=t(:,i);
                nodePos=n(:,nodeIndeces);
                e21local=e21El(nodePos);
                e21=max(e21,e21local);
            end
        end

        function e20=getE20Cons(obj)
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            %|u-Iu|_0<=E20 |u|_2, computed with Kobayashi formulae
            n=obj.nodes;
            t=obj.triangles;
            e20=0;
            for i=1:size(t,2)
                nodeIndeces=t(:,i);
                nodePos=n(:,nodeIndeces);
                e20local=e20El(nodePos);
                e20=max(e20,e20local);
            end
        end
        
        function setIsConvex(obj,convexity)
            % since the consistent convexity tester is sensitive to
            % rounding error, we also define the following method
            obj.isConvex=convexity;
        end
    
        function newMesh=uniformRefine(obj)
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            n=obj.nodes;
            t=obj.triangles;
            [newNodes,newTriangles]=refineMesh4(n,t);
            newMesh=P1Mesh(newNodes,newTriangles);
        end
    
        function visualize(obj,titleFig,path)
            n=obj.nodes;
            t=obj.triangles;
            t=[t;ones( 1,size(t,2) )];

            figure;
            pdemesh(n,[],t);

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
    
        function [K,C,R]=getEllipticMatrices(obj,dif,convection,reaction)
            %inputs:
            %   dif: 2*2 SPD matrix or positive real number, default:1
            %   convection: 2*1 vector, default: [0,0]'
            %   reaction: non-negative real number, default: 1

            %outputs:
            %   K_ij=\int (dif grad P_j)*grad P_i
            %   C_ij=\int (convection*grad P_j) P_i
            %   R_ij=\int reaction*P_j*P_i
            %   (including boundary nodes)
            arguments
                obj
                dif=1.0
                convection=[0.0,0.0]'
                reaction=1.0
            end
            % assemble K, M with built-in functions
            nodeTable=obj.nodes;
            triangleTable=obj.triangles;

            model=createpde();
            geometryFromMesh(model,nodeTable,triangleTable);
            mesh=model.Mesh;

            if isequal( size(dif),[1,1] )
                difInput=dif;
            end
            if isequal( size(dif),[2,2] )
                difInput=[dif(1,1);dif(1,2);dif(2,1);dif(2,2)];
            end
            % m u_tt + d u_t-div( c grad u)+au=f
            specifyCoefficients(model, 'm',reaction, 'd',0, ...
                'c',difInput, 'a',0, 'f',0);

            FEM = assembleFEMatrices(model,'KM');
            K=FEM.K;
            R=FEM.M;

            % triangle-wise assembly of C
            % assemble the convection matrix including boundary nodes on each
            % triangle
            numNodes = size(nodeTable,2);
            numTriangles = size(triangleTable,2);
        
            C=zeros(numNodes,numNodes);
            areas=zeros(numTriangles,1);
            for el=1:numTriangles
                areas(el)=area(mesh,el);
            end
            for el=1:numTriangles
                S=areas(el);
                % indeces of el
                n1=triangleTable(1,el);
                n2=triangleTable(2,el);
                n3=triangleTable(3,el);
                %coordinates of the verteces
                v1=nodeTable(:,n1);
                v2=nodeTable(:,n2);
                v3=nodeTable(:,n3);
        
                %gradients of the nodal basis of the verteces
                g1=grad(v1,v2,v3);
                g2=grad(v2,v1,v3);
                g3=grad(v3,v1,v2);
        
                % Compute contributions to the stiffness matrix C
                C(n2,n1) = C(n2,n1) + (g1' * convection) * S/3;
                C(n3,n1) = C(n3,n1) + (g1' * convection) * S/3;
                C(n1,n2) = C(n1,n2) + (g2' * convection) * S/3;
                C(n3,n2) = C(n3,n2) + (g2' * convection) * S/3;
                C(n1,n3) = C(n1,n3) + (g3' * convection) * S/3;
                C(n2,n3) = C(n2,n3) + (g3' * convection) * S/3;
            end
            C=sparse(C);
        end
    
        function diameter=getDiamater(obj)
            addpath(fullfile(fileparts(mfilename('fullpath')), '\Private'));
            nodeTable=obj.nodes;
            bNodes=obj.boundaryNodes;
            boundary=nodeTable(:,bNodes);
            diameter=computeDiameter(boundary);
        end
    end
end