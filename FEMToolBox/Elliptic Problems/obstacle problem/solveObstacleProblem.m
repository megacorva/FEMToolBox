function uh=solveObstacleProblem(initialGuess,dif,convection,reaction, ...
    force,obstacle,tolerance, useGPU)
    % finds u in K s.t. (Lu,v-u)>=(f,v-u),
    %  for all v in K={v in H_0^1:v>=obstacle}
    %  with the FEM scheme
    % finding u_h in K_h, s.t. (Lu_h, v_h-u_h)>=(f,v_h-u_h),
    % for all v_h in K_h={v in P1 FESpace:v>=b},
    % where 
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
    %   force: FEFunc
    %   obstacle:FEFunc
    %   tolerance: positive real number
    %   useGPU: boolean
    %----------------------------------------------------------------------
    % outputs:
    %   uh: FEFunc
    %----------------------------------------------------------------------
    % remark:
    %   the truncation error can be derived from Falk's lemma
    %   and regularity theories of variational inequalities,
    %   which involve a lot of variables that is not
    %   numerically computable,
    %   but we give a control on l^2 iterative error 
    %   in the coordinate space
    %----------------------------------------------------------------------

    projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
    addpath(genpath(projectRoot));

    mesh=initialGuess.mesh;
    internalNodes=mesh.internalNodes;
    % 1. get FE matrices---------------------------------------------------
    [K,C,R]=mesh.getEllipticMatrices(dif,convection,reaction);
    S=K+R;
    L=S+C;
    SI=S(internalNodes,internalNodes);
    LI=L(internalNodes,internalNodes);
    
    % 2. get the load vector-----------------------------------------------
    fNodalValues=force.nodalValues;
    [~,~,M]=mesh.getEllipticMatrices();
    F=M * fNodalValues;
    FI=F(internalNodes);

    % 3. interpolate the obstacle------------------------------------------
    obstacleNodal = obstacle.nodalValues;
    ob=obstacleNodal(internalNodes);

    % 4. iterative solution
    eigMin=eigs( SI ,1,'smallestreal');
    svMax=svds(LI,1,'largest');
    step=eigMin/svMax^2;
    contraction=sqrt(1-eigMin^2/svMax^2);
    % disp(contraction);
    eps=tolerance*(1-contraction)/contraction;

    % initial value
    u1=initialGuess.nodalValues(internalNodes);
    u2=u1+step*(FI-LI*u1);
    u2=clip(u2,ob,Inf);
    
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
                    u2=clip(u2,ob,Inf);
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
                u2=clip(u2,ob,Inf);
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