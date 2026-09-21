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
force=@(x,y)-10;
obstacle=@(x,y) -0.5+0.5*(1-x.^2).*(1-y.^2);

% setting the initial value------------------------------------------------
u0=zeros( size(mesh.nodes,2),1 );
initialGuess=FEFunc(mesh,u0);

% solve PDE and visualize the result and truncation error------------------

% solve on the initial coarse mesh
uh=solveObstacleProblem(initialGuess,dif,convection, ...
    reaction,force,obstacle,10^-6);

% % refine 1 time and solve on the refined mesh
% uh=uh.uniformRefine();
% 
% uh=solveObstacleProblem(uh,dif,convection, ...
%     reaction,force,obstacle,10^-8);

uh.visualize();