% ----------------------------------------------------------------------------
% LipidCA_main: Driver program for Lipid Cellular Automation
% Article: A. Gupta, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% DOI: https://doi.org/10.1039/C8SM02032A
% License: GNU GPL 3.0
% ----------------------------------------------------------------------------


clear all; clc;
%profile on

circleInitialRadius = 100;
circleMaxRadius = 150;
pinningProb = 0.02;  % Probability that cell becomes a pin
clusterProb = 0.0000; % Probability that pinned cell becomes a cluster pin
clusterThresh = 0.36; %log(0.5*circleInitialRadius/circleMaxRadius + 1);  % Probability that cells neighboring cluster cell join the cluster
percolationThresh = 0.7; % Probability that neighbor of a broken cluster cell becomes the root of a new offshoot cluster
bondStrength = 0.8; % pins fracture when total tension exceeds bondStrength
avalancheParameter = 20.0;  % Determines how much tension is released in offshoot cluster
tensionAdhesion = 0.2;  % tension at outer edge of lipid due to adhesion with substrate pulling radially outwards
tensionMLV = 0.1; % tension at center of membrane due to MLV opposing spread
tensionPin = 1.0; % tension parameter associated with pinning (usually 0 or 1)
clusterTensionOnOff = 1.0; % whether to include cluster tension (usually 0 or 1)
clusterTensionParameter = 0.75; % parameter inside cluster tension computation
radiusThresholdParameter = circleInitialRadius; % how many radial steps to look back when creating cluster chains.

numberOfIterations = circleMaxRadius-circleInitialRadius;
circleCurrentRadius = circleInitialRadius;
automationSize = 2*circleMaxRadius+1;
numCells = automationSize^2;


cellLocation = zeros(numCells, 3);
cellKey = zeros(automationSize, automationSize);
cell = 1;
for xLocation= -circleMaxRadius:circleMaxRadius
    col = xLocation+circleMaxRadius + 1;
    for yLocation = -circleMaxRadius:circleMaxRadius
        row = circleMaxRadius - yLocation + 1;
        cellKey(row, col) = cell;
        cellLocation(cell, 1) = xLocation;
        cellLocation(cell, 2) = yLocation;
        cellLocation(cell, 3) = sqrt(xLocation^2+yLocation^2);
        cell = cell + 1;
    end
end

cellStatus = zeros(numCells, 12);

for cell = 1:numCells
    cellStatus(cell, 6) = bondStrength;
end

cellNeighbors = getCellNeighbors(cellKey, automationSize, numCells);

clusterInfo = [];
numClusters = 0;
wettedArea = zeros(numberOfIterations,1);
fractureArea = zeros(numberOfIterations,1);

[cellStatus, clusterInfo, numClusters] = assignPinStatus(cellStatus, clusterInfo, cellLocation, circleCurrentRadius, numCells, numClusters, pinningProb, clusterProb);
[cellStatus, clusterInfo] = createClusterChain(cellStatus, clusterInfo, cellLocation, cellNeighbors, numClusters, clusterThresh, circleCurrentRadius, radiusThresholdParameter);
[wettedArea(1), fractureArea(1)] = computeAreas(cellStatus, numCells);

F(1) = getFigures(cellStatus, cellLocation, cellKey, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius);
for iteration = 1:numberOfIterations

    % Expand circle, assign new random numbers, finding pins and cluster pins
    circleCurrentRadius = circleCurrentRadius + 1;

    % Change cluster threshold
    %clusterThresh = log(0.5*circleCurrentRadius/circleMaxRadius + 1);

    [cellStatus, clusterInfo, numClusters] = assignPinStatus(cellStatus, clusterInfo, cellLocation, circleCurrentRadius, numCells, numClusters, pinningProb, clusterProb);

    % Invasion percolation loop over cluster pins
    [cellStatus, clusterInfo] = createClusterChain(cellStatus, clusterInfo, cellLocation, cellNeighbors, numClusters, clusterThresh, circleCurrentRadius, radiusThresholdParameter);

    % Compute tensions
    cellStatus = computeTension(cellStatus, cellLocation, circleCurrentRadius, numCells, circleMaxRadius, tensionAdhesion, tensionMLV, tensionPin);

    % Modify bond strength
    [cellStatus, clusterInfo] = modifyClusterTension(cellStatus, clusterInfo, numClusters, bondStrength, clusterTensionOnOff, clusterTensionParameter);

    % Check for fracture
    [cellStatus, clusterInfo, numClusters, F] = checkFracture(F, cellStatus, clusterInfo, cellKey, cellLocation, cellNeighbors, numClusters, numCells, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius, avalancheParameter);

    % Compute areas
    [wettedArea(iteration), fractureArea(iteration)] = computeAreas(cellStatus, numCells);

    % Get figure
    F(end+1) = getFigures(cellStatus, cellLocation, cellKey, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius);
end
