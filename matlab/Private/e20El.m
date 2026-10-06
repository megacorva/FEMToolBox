function err=e20El(el)
    %el=[x1,x2,x3;
    % y1;y2,y3]

    %Computing edge length^2
    L12=(el(1,1)-el(1,2))^2+(el(2,1)-el(2,2))^2;
    L13=(el(1,1)-el(1,3))^2+(el(2,1)-el(2,3))^2;
    L23=(el(1,2)-el(1,3))^2+(el(2,2)-el(2,3))^2;

    area=polyarea(el(1,:),el(2,:));


    Kobayashi20=(L12*L13+L12*L23+L13*L23)/83-...
        (1/24)*((L12*L23*L13)/(L12+L13+L23)+area^2);
    err=sqrt(Kobayashi20);
    
end