% solves an obstacle problem, and visualize the result

projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(projectRoot));

% generating a Delauney mesh on the unit square----------------------------
partition=20;
[nodes, ~, triangles] = poimesh(@squareg, partition, partition);
triangles=triangles(1:3,:);
mesh=P1Mesh(nodes,triangles);

% setting parameters of the obstacle problem-------------------------------
dif=1;
convection=[0;0];
reaction=0;
force=@(x,y) -0.4+0.5*(1-x.^2).*(1-y.^2);

truncError=getZeroObstacleTruncErr(mesh,dif,convection,reaction,force);
consistency=truncError*2;
fprintf(['the consistency error from the problem with interpolated ' ...
    'force is at most %.3e'],consistency);

% setting the initial value------------------------------------------------
initialGuess=solveEllipticPDE(mesh,dif,convection,reaction,force);

% solve and visualize the result and truncation error----------------------

uh=solveZeroObstacleProblem(initialGuess,dif,convection, ...
    reaction,force,truncError);

uh.visualize();