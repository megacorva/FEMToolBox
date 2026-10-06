function diameter=computeDiameter(boundary)
    % inputs:
    %[x1,x2,...,x_n;
    %y1,y2,...,y_n]
    % sorted in the counterclockwise direction
    diameter=0.0;
    numPoints = size(boundary, 2);
    for i = 1:numPoints-1
        for j = i+1:numPoints
            dist = norm(boundary(:, i) - boundary(:, j));
            diameter = max(diameter, dist);
        end
    end
end