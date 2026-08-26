function isConvex=testConvexity(boundary)
    % inputs:
    %[x1,x2,...,x_n;
    %y1,y2,...,y_n]
    % sorted in the counterclockwise direction
    isConvex=true;
    numNodes=size(boundary,2);
    for i=1:numNodes
        previous=i-1;
        current=i;
        next=i+1;
        if i==1
            previous=numNodes;
        end
        if i==numNodes
            next=1;
        end
        %the vector from the previous node to the current node
        vecPC=boundary(:,current)-boundary(:,previous);
        %the vector from the current node to the next node
        vecPN=boundary(:,next)-boundary(:,previous);
        
        % Compute the cross product to determine the turn direction
        crossProd = vecPC(1) * vecPN(2) - vecPC(2) * vecPN(1);
        if crossProd < 0
            isConvex = false; % Not convex if the cross product is negative
            return; % Exit the function early
        end
    end
end