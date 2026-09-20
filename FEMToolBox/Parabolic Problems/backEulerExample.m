% solves a parabolic PDE with backward Euler method, 
% and visualize the result

addpath(fullfile(fileparts(mfilename('fullpath')), '\..'));

% generating a Delauney mesh on the unit square----------------------------
partition=40;
[nodes, ~, triangles] = poimesh(@squareg, partition, partition);
triangles=triangles(1:3,:);
mesh=P1Mesh(nodes,triangles);

% Define time step and total time------------------------------------------
dt = 0.01; 
totalTime = 1.0; 
nTime = totalTime / dt;

% setting parameters of the PDE--------------------------------------------
dif=1;
convection=[0;0];
reaction=0;
initial=@(x,y)0;
force=@(x,y)(1-x.^2).*(1-y.^2);

% solve and visualize------------------------------------------------------
uh=backEuler(mesh,initial,dif,convection,reaction,force,dt,nTime);
uh.visualize(20,'Heat Equation','heatEquation.mp4');

