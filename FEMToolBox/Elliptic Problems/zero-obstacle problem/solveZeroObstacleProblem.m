function uh=solveZeroObstacleProblem(initialGuess,dif,convection,reaction, ...
    force,H1tolerance, useGPU)
    % finds u in K s.t. (Lu,v-u)>=(f_h,v-u),
    %  for all v in K={v in H_0^1:v>=0}
    %  with the FEM scheme
    % finding u_h in K_h, s.t. (Lu_h, v_h-u_h)>=(f_h,v_h-u_h),
    % for all v_h in K_h={v in P1 FESpace:v>=0},
    % where 
    % Lu=-div( dif*grad u) + convection*grad u + reaction* u,
    % iterative error bound:
    % |uh-uh*|_1<=tolerance,
    % where uh* is the exact solution of the FE variational inequality.
    %----------------------------------------------------------------------
    % inputs:
    %   initialGuess: FEFunc
    %   dif: 2*2 SPD matrix, or positive real number
    %   convection: 2*1 vector
    %   reaction: non-negative real number
    %   force: FEFunc
    %   H1tolerance: positive real number
    %   useGPU: boolean
    %----------------------------------------------------------------------
    % outputs:
    %   uh: FEFunc
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
    [~,~,M]=mesh.getEllipticMatrices();
    fNodalValues=force.nodalValues;
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
    
    if ~isequal(u1,u2)
        % GPU parallelization
        if useGPU && canUseGPU
            LI=gpuArray(LI);
            FI=gpuArray(FI);
            u1=gpuArray(u1);
            u2=gpuArray(u2);
            
            itNum=0;
            while norm(u1-u2)>eps
                for i=1:50
                    u1=u2;
                    u2=u2+step*(FI-LI*u2);
                    u2=clip(u2,0,Inf);
                end
                itNum=itNum+50;
            end
            fprintf('%d GPU iterations till convergence.\n' ...
            ,itNum);
            u2=gather(u2);
        else
            itNum=0;
            while norm(u1-u2)>=eps
                u1=u2;
                u2=u2+step*(FI-LI*u2);
                u2=clip(u2,0,Inf);
                itNum=itNum+1;
            end
            fprintf('%d CPU iterations till convergence.\n', ...
            itNum);
        end
    end
    
    u=zeros( size(mesh.nodes,2),1 );
    u(internalNodes)=u2;
    uh=FEFunc(mesh,u);
end