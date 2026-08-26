%find the gradient of the nodal basis at des in the triangle
%(des,v1,v2)
function g=grad(des,v1,v2)
    t1=v1-des;
    t2=v2-des;
    g = (1/det([t1,t2]))*[t1(2)-t2(2);t2(1)-t1(1)];
end