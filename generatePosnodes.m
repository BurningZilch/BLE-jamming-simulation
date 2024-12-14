function Posnodes = generatePosnodes(pointDistance, xRange, yRange, zRange)
    % Generates evenly spread positions for Posnodes in a 3D area based on
    % specified point distance.
    % Inputs:
    % - pointDistance: Desired distance between each node
    % - xRange: [minX, maxX] range for x-coordinates
    % - yRange: [minY, maxY] range for y-coordinates
    % - zRange: [minZ, maxZ] range for z-coordinates

    % Calculate the number of nodes per dimension based on point distance
    numNodesX = floor((xRange(2) - xRange(1)) / pointDistance) + 1;
    numNodesY = floor((yRange(2) - yRange(1)) / pointDistance) + 1;
    numNodesZ = floor((zRange(2) - zRange(1)) / pointDistance) + 1;
    
    % Generate grid points within the specified ranges
    xPoints = linspace(xRange(1), xRange(2), numNodesX);
    yPoints = linspace(yRange(1), yRange(2), numNodesY);
    zPoints = linspace(zRange(1), zRange(2), numNodesZ);

    % Initialize coordinates for nodes
    totalNodes = numNodesX * numNodesY * numNodesZ;
    Posnodes = zeros(3, totalNodes);

    % Set a random offset range (e.g., 10% of the distance to spread points evenly)
    offsetFactor = 0.1;

    % Fill Posnodes with coordinates from the 3D grid, adding random offsets
    nodeIndex = 1;
    for i = 1:numNodesX
        for j = 1:numNodesY
            for k = 1:numNodesZ
                % Add random offsets to x, y, and z coordinates
                randomOffsetX = pointDistance * offsetFactor * (rand - 0.5) * 2;
                randomOffsetY = pointDistance * offsetFactor * (rand - 0.5) * 2;
                randomOffsetZ = pointDistance * offsetFactor * (rand - 0.5) * 2;

                % Assign positions with offsets
                Posnodes(:, nodeIndex) = [xPoints(i) + randomOffsetX; 
                                          yPoints(j) + randomOffsetY; 
                                          zPoints(k) + randomOffsetZ];
                nodeIndex = nodeIndex + 1;
            end
        end
    end
end
