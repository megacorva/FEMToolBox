function H1Err=getZeroObstacleTruncErr(mesh,dif,convection,reaction, ...
    force)
    % gets an H1 error bound in
    % finding u in K s.t. (Lu,v-u)>=(f_h,v-u),
    %  for all v in K={v in H_0^1:v>=0}
    %  with the FEM scheme
    % finding u_h in K_h, s.t. (Lu_h, v_h-u_h)>=(f_h,v_h-u_h),
    % for all v_h in K_h={v in P1 FESpace:v>=0},
    % where 
    % Lu=-div( dif*grad u) + convection*grad u + reaction* u,
    % iterative error bound:
    % |uh-uh*|<=tolerance,
    % where uh* is the exact solution of the FE variational inequality.
    %----------------------------------------------------------------------
    % inputs:
    %   mesh:P1Mesh
    %   dif: 2*2 SPD matrix, or positive real number
    %   convection: 2*1 vector
    %   reaction: non-negative real number
    %   force: FEFunc
    %----------------------------------------------------------------------
    % outputs:
    %   H1Err: positive real number
    %----------------------------------------------------------------------

    H1Err=NaN;
    % 1. getting the eigenvalues of dif and Poincare constant--------------
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
    
    % 2. getting the coercity constant and operator bound------------------
    opCoercity=eigDifMin + min(0,reaction*CP^2);
    
    if opCoercity>0

        % 3. getting regularity constant and operator bound----------------
        regularity=( 1 + ( norm(convection) + abs(reaction)*CP ) ...
                    *CP/opCoercity )/eigDifMin;
        opBound=eigDifMax + norm(convection)*CP + abs(reaction*CP^2);

        % 4. getting interpolation error constants-------------------------
        I0=mesh.getE20Cons();
        I1=mesh.getE21Cons();
        
        % 5. getting L2 norm of f_h----------------------------------------
        [~,~,M]=mesh.getEllipticMatrices();
        fNodalValues=force.nodalValues;
        F0=sqrt(fNodalValues'*M*fNodalValues);

        H1Err=F0*( I0*2*regularity/opCoercity + ...
           ( I1 * opBound *regularity / opCoercity )^2 );
    end
end