function err=e21El(el)
    %el=[x1,x2,x3;
    % y1;y2,y3]

    %Computing edge length^2
    L1Sq=(el(1,1)-el(1,2))^2+(el(2,1)-el(2,2))^2;
    L2Sq=(el(1,1)-el(1,3))^2+(el(2,1)-el(2,3))^2;
    L3Sq=(el(1,2)-el(1,3))^2+(el(2,2)-el(2,3))^2;

    area=polyarea(el(1,:),el(2,:));


    Kobayashi2=L1Sq*L2Sq*L3Sq/(16*area^2) ...
    -(L1Sq+L2Sq+L3Sq)/30 ...
    -(area^2/5)*(1/L1Sq+1/L2Sq+1/L3Sq);
    err=sqrt(Kobayashi2);
    
end