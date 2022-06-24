% ----------------------------------------------------------------------------
% clusterAvalanche: Fracture routine for cluster pins
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function  [cellStatus, clusterInfo, numClusters] = clusterAvalanche(F, cellStatus, clusterInfo, cellKey, cellLocation, cellNeighbors, numClusters, numCells, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius, avalancheParameter, radiusThresh, clusterBreakingTension, cluster, numCellsInCluster)

randSearch = 1;
cellSearch = clusterInfo(4,cluster);
for clusterCell = 1:numCellsInCluster
    cell = clusterInfo(clusterCell+3, cluster);
    cellStatus(cell,10) = 1; % break all pins in chain
end
for clusterCell = 1:numCellsInCluster
    cell = clusterInfo(clusterCell+3, cluster);
    % loop over neighbors, keep tally of lowest R1
    for neighbor = 1:8
        if cellNeighbors(cell, neighbor) ~= -1
            if (cellStatus(cellNeighbors(cell, neighbor), 1) == 1)&&(cellStatus(cellNeighbors(cell, neighbor), 2) < randSearch)&&(cellStatus(cellNeighbors(cell, neighbor), 4) == 0)
                randSearch = cellStatus(cellNeighbors(cell, neighbor), 2);
                cellSearch = cellNeighbors(cell, neighbor);
            end
            if (cellStatus(cellNeighbors(cell, neighbor), 5) == 1)&&(cellStatus(cellNeighbors(cell, neighbor), 10) == 0) % if neighbor is an unbroken cluster cell
                neighborClusterNumber = cellStatus(cellNeighbors(cell, neighbor), 12);
                numCellsInNeighborCluster = clusterInfo(1, neighborClusterNumber);
                if (neighborClusterNumber ~= cluster)
                    [cellStatus, clusterInfo, numClusters] = clusterAvalanche(F, cellStatus, clusterInfo, cellKey, cellLocation, cellNeighbors, numClusters, numCells, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius, avalancheParameter, radiusThresh, clusterBreakingTension, neighborClusterNumber, numCellsInNeighborCluster);
                end
            end
        end
    end
end
if (randSearch < percolationThresh)&&(cellStatus(cellSearch, 5) == 0) %create new offshoot cluster
    cellStatus(cellSearch, 4) = 1;
    cellStatus(cellSearch, 5) = 1;
    numClusters = numClusters+1;
    cellStatus(cellSearch, 12) = numClusters;
    clusterInfo(1, numClusters) = 1; % numCells in cluster
    clusterInfo(2, numClusters) = 0; % cluster is initially unbroken
    clusterInfo(3, numClusters) = 1; % cluster is an offshoot
    clusterInfo(4, numClusters) = cellSearch; % cell number
    [cellStatus, clusterInfo] = findClusterPinNeighbors(cellStatus, clusterInfo, cellLocation, cellNeighbors, numClusters, cellSearch, radiusThresh, clusterThresh);
    for clusterCell = 1:clusterInfo(1, numClusters) % assign Treleased to all offshoot cluster pins
        offshootCell = clusterInfo(clusterCell+3, numClusters);
        cellStatus(offshootCell, 9) = -(clusterBreakingTension - cellStatus(offshootCell, 6))*avalancheParameter;
    end

end
