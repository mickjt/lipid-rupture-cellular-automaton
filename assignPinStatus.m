% ----------------------------------------------------------------------------
% assignPinStatus: Determine whether cells are pins or not
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function [cellStatus, clusterInfo, numClusters] = assignPinStatus(cellStatus, clusterInfo, cellLocation, circleCurrentRadius, numCells, numClusters, pinningProb, clusterProb)

for cell = 1:numCells
    if (cellLocation(cell, 3) <= circleCurrentRadius)&&(cellStatus(cell,1) ~= 1) % cell within the lipid and not already on
        cellStatus(cell, 1) = 1;
        cellStatus(cell, 2) = rand;
        cellStatus(cell, 3) = rand;
        if cellStatus(cell,2) < pinningProb
            cellStatus(cell, 4) = 1;
            if cellStatus(cell, 3) < clusterProb
                cellStatus(cell, 5) = 1;
                numClusters = numClusters+1;
                cellStatus(cell, 12) = numClusters;
                clusterInfo(1, numClusters) = 1; % numCells in cluster
                clusterInfo(2, numClusters) = 0; % cluster is initially unbroken
                clusterInfo(3, numClusters) = 0; % cluster is not an offshoot
                clusterInfo(4, numClusters) = cell; % cell number
            end
        end
    end
end
