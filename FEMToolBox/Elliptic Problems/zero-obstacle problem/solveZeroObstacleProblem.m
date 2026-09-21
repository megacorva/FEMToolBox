function uh=solveZeroObstacleProblem(initialGuess,dif,convection,reaction, ...
    force,H1tolerance)
    % finds u in K s.t. (Lu,v-u)>=(f_h,v-u),
    %  for all v in K={v in H_0^1:v>=0}
    %  with the FEM scheme
    % finding u_h in K_h, s.t. (Lu_h, v_h-u_h)>=(f_h,v_h-u_h),
    % for all v_h in K_h={v in P1 FESpace:v>=0},
    % where 
    % f_h is the Lagrange interpolation of force,
    % Lu=-div( dif*grad u) + convection*grad u + reaction* u,
    % iterative error bound:
    % |uh-uh*|<=tolerance,
    % where uh* is the exact solution of the FE variational inequality.
    %----------------------------------------------------------------------
    % inputs:
    %   initialGuess: FEFunc
    %   dif: 2*2 SPD matrix, or positive real number
    %   convection: 2*1 vector
    %   reaction: non-negative real number
    %   force=@(x,y)...
    %   H1tolerance: positive real number
    %----------------------------------------------------------------------
    % outputs:
    %   ub: FEFunc
    %----------------------------------------------------------------------


    projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
    addpath(genpath(projectRoot));

    mesh=initialGuess.mesh;
    internalNodes=mesh.internalNodes;
    % 1. get FE matrices---------------------------------------------------
    [H,~,~]=mesh.getEllipticMatrices(1,[0;0],0);
    [K,C,R]=mesh.getEllipticMatrices(dif,convection,reaction);
    S=K+R;
    L=S+C;
    HI=H(internalNodes,internalNodes);
    SI=S(internalNodes,internalNodes);
    LI=L(internalNodes,internalNodes);
    
    % 2. get the load vector-----------------------------------------------
    xNodes=mesh.nodes(1,:);
    yNodes=mesh.nodes(2,:);
    [~,~,M]=mesh.getEllipticMatrices();
    fNodalValues=arrayfun(force, xNodes, yNodes)';
    F=M * fNodalValues;
    FI=F(internalNodes);

    % 3. iterative solution
    eigHMax=eigs ( HI, 1, 'largestreal' );
    eigMin=eigs( SI ,1,'smallestreal');
    svMax=svds(LI,1,'largest');
    step=eigMin/svMax^2;
    contraction=sqrt(1-eigMin^2/svMax^2);
    % disp(contraction);
    eps=H1tolerance*(1-contraction)/( contraction * eigHMax );

    % initial value
    u1=initialGuess.nodalValues(internalNodes);
    u2=u1+step*(FI-LI*u1);
    u2=clip(u2,0,Inf);

    % parallelization
    if canUseGPU
        LI=gpuArray(LI);
        FI=gpuArray(FI);
        u1=gpuArray(LI);
        u2=gpuArray(LI);
    end

    while norm(u1-u2)>=eps
        u1=u2;
        u2=u2+step*(FI-LI*u2);
        u2=clip(u2,0,Inf);
    end
    
    u=zeros( size(mesh.nodes,2),1 );

    if canUseGPU
        u2=gather(u2);
    end

    u(internalNodes)=u2;
    uh=FEFunc(mesh,u);
end