function [uh,H1Error,H0Error]=solveEllipticPDE(mesh,dif,convection,reaction,f)
    % solves Lu=f_h with FEM, 
    % where f_h is the Lagrange interpolation of f,
    % Lu=-div( dif*grad u) + convection*grad u + reaction* u,
    % u=0 on the boundary of the mesh
    %----------------------------------------------------------------------
    %inputs:
    %   mesh: P1Mesh
    %   dif: 2*2 SPD matrix, or positive real number
    %   convection: 2*1 vector
    %   reaction: non-negative real number
    %   f=@(x,y)...
    %----------------------------------------------------------------------
    %outputs:
    %   uh: FEFunc
    %   H1Error: |uh-u|_1
    %   H0Error: ||uh-u||_0
    %   where Lu=f_h, f_h is the Lagrange interpolation of f onto the mesh
    %----------------------------------------------------------------------
    
    addpath(fullfile(fileparts(mfilename('fullpath')), '\..'));
    internalNodes=mesh.internalNodes;
    % 1. get FE matrices---------------------------------------------------
    [K,C,R]=mesh.getEllipticMatrices(dif,convection,reaction);
    L=K+C+R;
    LI=L(internalNodes,internalNodes);
    
    % 2. get the load vector-----------------------------------------------
    xNodes=mesh.nodes(1,:);
    yNodes=mesh.nodes(2,:);
    [~,~,M]=mesh.getEllipticMatrices();
    fNodalValues=arrayfun(f, xNodes, yNodes)';
    F=M * fNodalValues;
    FI=F(internalNodes);

    % 3. find the nodal values of the FE solution--------------------------
    uI=LI\FI;

    % 4. assemble uh-------------------------------------------------------
    u=zeros( size(mesh.nodes,2),1 );
    u(internalNodes)=uI;
    uh=FEFunc(mesh,u);

    % 5. estimate H1Error and H0Error--------------------------------------
    H1Error=NaN;
    H0Error=NaN;
    
    % get an estimate of truncation error by Cea's lemma, duality argument
    % and Miranda-Talenti inequality
    if mesh.isConvex
        if isequal( size(dif),[1,1] )
            eigDifMin=dif;
            eigDifMax=dif;
        else
            eigsDif=eig(dif);
            eigDifMin = min(eigsDif);
            eigDifMax = max(eigsDif);
        end

        % a simple upper bound of the Poincare constant
        % |u|_0<= CP|u|_1
        CP=mesh.getDiamater()/ ( sqrt(2)*pi );

        opCoercity=eigDifMin + min(0,reaction*CP^2);
        if opCoercity>0
            opBound=eigDifMax + norm(convection)*CP + max(0,reaction*CP^2);

            % L^2 norm of f_h
            f0=sqrt(fNodalValues'*F);
            
            % upper bound of H^1 seminorm of u
            u1=f0*CP/opCoercity;

            % estimate |u|_2 by Miranda-Talenti
            u2=(  f0 + ( norm(convection) + abs(reaction)*CP )*u1  )/eigDifMin;

            % interpolation error constants
            e21=mesh.getE21Cons();

            % estimate |u-uh|_1 by Cea's lemma
            H1Error=(opBound/opCoercity) * e21* u2;

            % estimate ||u-uh||_0 using duality argument
            H0Error=opBound*e21*H1Error/eigDifMin;
        end
    end
end