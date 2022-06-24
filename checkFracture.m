% ----------------------------------------------------------------------------
% checkFracture: Determine whether pinned cells have broken and call appropriate subroutines
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function  [cellStatus, clusterInfo, numClusters, F] = checkFracture(F, cellStatus, clusterInfo, cellKey, cellLocation, cellNeighbors, numClusters, numCells, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius, avalancheParameter)

radiusThresh = 0;
for cluster = 1:numClusters
    numCellsInCluster = clusterInfo(1, cluster);
    clusterBreakingTension = 0;
    if clusterInfo(2, cluster) == 0 % if cluster is not broken
        % Check to see if cluster breaks
        for clusterCell = 1:numCellsInCluster
            cell = clusterInfo(clusterCell+3, cluster);
            totalTension = cellStatus(cell, 7) + cellStatus(cell, 8) + cellStatus(cell, 9);
            if (totalTension > cellStatus(cell,6)) % if tension > bondstrength
                clusterBreakingTension = totalTension;
                cellStatus(cell,10) = 1;
                clusterInfo(2, cluster) = 1;
                break;
            end
        end

        % if cluster breaks
        if clusterInfo(2, cluster) == 1
            [cellStatus, clusterInfo, numClusters] = clusterAvalanche(F, cellStatus, clusterInfo, cellKey, cellLocation, cellNeighbors, numClusters, numCells, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius, avalancheParameter, radiusThresh, clusterBreakingTension, cluster, numCellsInCluster);
        end
    end
end


test = 1:numCells;
test_rand = test(randperm(length(test)));
for cell = test_rand %1:numCells %make this randomly ordered
    % If it's an unbroken pin: break a one radius circle around it, check
    % for pins, break those new pins with same radius
    if (cellStatus(cell,4) == 1) && (cellStatus(cell,5) == 0)  % If cell is a pin but not a cluster pin
        if (cellStatus(cell,7) > cellStatus(cell,6)) || (cellStatus(cell,10) == 1) % if tension > bondStrength or pin is already broken
            cellStatus(cell,10) = 1;
            [F, cellStatus] = breakNeighbors(F, cellStatus, cellKey, cellLocation, numCells, automationSize, cell, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius);
        end
    end

end
