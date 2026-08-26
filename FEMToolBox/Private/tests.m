% Define the mesh parameters (nx and ny are number of divisions)
nx = 2; ny = 2;
[p, e, t] = poimesh(@squareg, nx, ny);

p(:,6)=[0.5;0];
t=t(1:3,:);

bNodes=findBoundaryNodes(p,t);
boundary=p(:,bNodes);

disp(testConvexity(boundary));

[p,t]=refineMesh4(p,t);

model=createpde();
geometryFromMesh(model,p,t);
figure
pdemesh(model);