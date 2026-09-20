% solves an elliptic PDE, visualize the result together with explicit
% truncation error bounds

addpath(fullfile(fileparts(mfilename('fullpath')), '\..'));

% generating a Delauney mesh on the unit square----------------------------
partition=100;
[nodes, ~, triangles] = poimesh(@squareg, partition, partition);
triangles=triangles(1:3,:);
mesh=P1Mesh(nodes,triangles);

% setting parameters of the PDE--------------------------------------------
dif=[1,0;0,0.5];
convection=[2;1];
reaction=0;
f=@(x,y)1;

% solve PDE and visualize the result and truncation error------------------
[uh,H1Error,H0Error]=solveEllipticPDE(mesh,dif,convection,reaction,f);

titleFig= "H1 error bound: "+H1Error +", L2 error bound:"+H0Error;

uh.visualize(titleFig);