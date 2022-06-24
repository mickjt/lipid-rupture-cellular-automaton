% ----------------------------------------------------------------------------
% createClusterChain: Create a chain of cluster pin sites
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function [cellStatus, clusterInfo] = createClusterChain(cellStatus, clusterInfo, cellLocation, cellNeighbors, numClusters, clusterThresh, circleCurrentRadius, radiusThresholdParameter)

radiusThresh = circleCurrentRadius - radiusThresholdParameter; % look back only radiusThresholdParameter radial steps
for cluster = 1:numClusters
    numCellsInCluster = clusterInfo(1, cluster);
    for clusterCell = 1:numCellsInCluster
        cell = clusterInfo(clusterCell+3, cluster);
        [cellStatus, clusterInfo] = findClusterPinNeighbors(cellStatus, clusterInfo, cellLocation, cellNeighbors, cluster, cell, radiusThresh, clusterThresh);
    end
end
