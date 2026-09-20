% force=@(t,x,y)...
% initial=@(x,y)...
function uh=backEuler(mesh, initial,dif,convection,reaction, ...
    force, dt, nTime)
    % solves u_t+Lu=f_h, u(0)=initial with FEM backward Euler method,
    % where f_h is the Lagrange interpolation of f,
    % Lu=-div( dif*grad u) + convection*grad u + reaction* u,
    % u=0 on the boundary of the mesh
    %----------------------------------------------------------------------
    % inputs:
        % mesh:P1Mesh
        % initial=@(x,y)...Initial value of solution
        % dif: 2*2 SPD matrix or positive scalar
        % convection:2*1 vector
        % reaction: scalar large enough so that L is coercive
        % force=@(t,x,y)...force
        % dt: positive scalar,time step length
        % nTime:integer, number of time samples
    %----------------------------------------------------------------------

    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
    
    internalNodes=mesh.internalNodes;

    % 1. initialize u and samplingTime, u=[u(0),...,u(tN)]-----------------
    u=zeros( size(mesh.nodes,2),nTime+1 );
    samplingTime=0:dt:nTime*dt;
    
    % Set initial conditions for internal nodes
    XI=mesh.nodes(1,:);
    YI=mesh.nodes(2,:);
    u(:, 1) = arrayfun(initial,XI,YI)'; 

    % 2. compute FE matrices-----------------------------------------------
    [K,C,R]=mesh.getEllipticMatrices(dif,convection,reaction);
    L=K+C+R;
    LI=L(internalNodes,internalNodes);
    [~,~,M]=mesh.getEllipticMatrices();
    MI=M(internalNodes,internalNodes);
    E=MI+dt*LI;
    [LL,UL]=lu(E);

    % 3. compute load vector,f=[f(t1),...f(tN)]----------------------------
    f = zeros( length(internalNodes), nTime ); 
    for n = 1:nTime
        % Compute load vector at each time step
        fNow=@(x,y)f(n*dt,x,y);
        fNodalValues = arrayfun(force, XI, YI)'; 
        load=M*fNodalValues;
        f(:,n)=load(internalNodes);
    end
    % time step
    for n=1:nTime
        now=u(internalNodes,n);
        w=MI*now+dt*f(:,n);
        v=UL\(LL\w);
        u(internalNodes,n+1)=v;
    end
    uh=TFEFunc(mesh,samplingTime,u);
end