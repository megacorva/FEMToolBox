function uh=solveObstacleProblem(initialValue,dif,convection,reaction, ...
    force,obstacle,tolerance)
    % finds u in K s.t. (Lu,v-u)>=(f_h,v-u),
    %  for all v in K={v in H_0^1:v>=obstacle}
    %  with the FEM scheme
    % finding u_h in K_h, s.t. (Lu_h, v_h-u_h)>=(f_h,v_h-u_h),
    % for all v_h in K_h={v in P1 FESpace:v>=b_h},
    % where b_h is the Lagrange interpolation of obstacle,
    % f_h is the Lagrange interpolation of force,
    % Lu=-div( dif*grad u) + convection*grad u + reaction* u,
    % iterative error bound:
    % |uh-uh*|<=tolerance,
    % where uh* is the exact solution of the FE variational inequality.
    %----------------------------------------------------------------------
    % inputs:
    %   initialValue: FEFunc
    %   dif: 2*2 SPD matrix, or positive real number
    %   convection: 2*1 vector
    %   reaction: non-negative real number
    %   force=@(x,y)...
    %   obstacle=@(x,y)...
    %   tolerance: positive real number
    %----------------------------------------------------------------------
    % outputs:
    %   ub: FEFunc
    %----------------------------------------------------------------------
    % remark:
    %   the truncation error can be derived from Falk's lemma, 
    %   since |obstacle|_2 is not generally numerically computable,
    %   we do not return the truncation error,
    %   but we give a control on l^2 iterative error 
    %   in the coordinate space
    %----------------------------------------------------------------------
    mesh=initialValue.mesh;
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
    fNodalValues=arrayfun(force, xNodes, yNodes)';
    F=M * fNodalValues;
    FI=F(internalNodes);

    % 3. interpolate the obstacle------------------------------------------
    obstacleFull = arrayfun(obstacle, xNodes, yNodes)';
    ob=obstacleFull(internalNodes);

    % 4. iterative solution
    svMin=svds(LI,1,'smallest');
    svMax=svds(LI,1,'largest');
    step=svMin/svMax^2;
    contraction=sqrt(1-svMin^2/svMax^2);
    % disp(contraction);
    eps=tolerance*(1-contraction)/contraction;

    % initial value
    u1=initialValue.nodalValues(internalNodes);
    u2=u1+step*(FI-LI*u1);
    u2=clip(u2,ob,Inf);

    while norm(u1-u2)>=eps
        u1=u2;
        u2=u2+step*(FI-LI*u2);
        u2=clip(u2,ob,Inf);
    end
    u=zeros( size(mesh.nodes,2),1 );
    u(internalNodes)=u2;
    uh=FEFunc(mesh,u);
end