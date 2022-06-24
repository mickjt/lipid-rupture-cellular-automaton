% ----------------------------------------------------------------------------
% findClusterPinNeighbors: Recursively determine which neighbor cells should join a given cluster
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function [cellStatus, clusterInfo] = findClusterPinNeighbors(cellStatus, clusterInfo, cellLocation, cellNeighbors, cluster, cell, radiusThreshold, clusterThresh)

for neighbor = 1:8
    if cellNeighbors(cell, neighbor) ~= -1
        if (cellStatus(cellNeighbors(cell, neighbor), 1) == 1)&&(cellStatus(cellNeighbors(cell, neighbor), 2) < clusterThresh) && (cellStatus(cellNeighbors(cell, neighbor), 5) ~= 1)&&(cellLocation(cell,3) > radiusThreshold)
            cellStatus(cellNeighbors(cell, neighbor), 4) = 1;
            cellStatus(cellNeighbors(cell, neighbor), 5) = 1;
            cellStatus(cellNeighbors(cell, neighbor), 12) = cluster;
            clusterInfo(1, cluster) = clusterInfo(1, cluster) + 1;
            clusterInfo(clusterInfo(1, cluster)+3, cluster) = cellNeighbors(cell, neighbor);
            [cellStatus, clusterInfo]  = findClusterPinNeighbors(cellStatus, clusterInfo, cellLocation, cellNeighbors, cluster, cellNeighbors(cell, neighbor), radiusThreshold, clusterThresh);
        end
    end
end
