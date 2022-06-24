% ----------------------------------------------------------------------------
% modifyClusterTension: reduce cluster tension based on size and CTP
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function [cellStatus, clusterInfo] = modifyClusterTension(cellStatus, clusterInfo, numClusters, bondStrength, clusterTensionOnOff, clusterTensionParameter)

for cluster = 1:numClusters
    numCellsInCluster = clusterInfo(1, cluster);
    for clusterCell = 1:numCellsInCluster
        cell = clusterInfo(clusterCell+3, cluster);
        cellStatus(cell,8) = -(numCellsInCluster*(clusterTensionParameter)^(numCellsInCluster-1) - 1) * bondStrength * clusterTensionOnOff;
    end
end
