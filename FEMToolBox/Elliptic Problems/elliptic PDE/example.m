% solves an elliptic PDE, visualize the result together with explicit
% truncation error bounds

projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(projectRoot));

% getting the hexagon mesh-------------------------------------------------
data = load(fullfile('FEMToolBox', 'Meshes', 'hexagon.mat'));
mesh = P1Mesh(data.nodes, data.triangles);
mesh=mesh.uniformRefine();

% setting parameters of the PDE--------------------------------------------
dif=[1,0;0,0.5];
convection=[0.5;0];
reaction=0;
fExact=@(x,y)1;

% solve PDE and visualize the result and truncation error------------------
fInterpolated=FEInterpolate(mesh,fExact);
[uh,H1Error,H0Error]=solveEllipticPDE(mesh,dif,convection,reaction, ...
    fInterpolated);

titleFig= "H1 error bound: "+H1Error +", L2 error bound:"+H0Error;

uh.visualize(titleFig);