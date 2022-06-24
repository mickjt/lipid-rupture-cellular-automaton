% ----------------------------------------------------------------------------
% breakNeighbors: Fracture routine for dilute pins
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function [F, cellStatus] = breakNeighbors(F, cellStatus, cellKey, cellLocation, numCells, automationSize, pinCell, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius)

cellStatus(pinCell, 11) = cellStatus(pinCell, 11) + 1;   % increment fracture radius

for cell = 1:numCells
    distanceFromPin = sqrt( (cellLocation(cell,1)-cellLocation(pinCell,1))^2 + (cellLocation(cell,2)-cellLocation(pinCell,2))^2 );
    if (distanceFromPin <= cellStatus(pinCell, 11)) && (cellStatus(cell, 10) == 0)
        cellStatus(cell, 10) = 1;
        if cellStatus(cell, 4) == 1 % cell is a pin
            for i = 1:cellStatus(pinCell, 11)
                 [F, cellStatus] = breakNeighbors(F, cellStatus, cellKey, cellLocation, numCells, automationSize, cell, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius);
                 F(end+1) = getFigures(cellStatus, cellLocation, cellKey, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius);
            end
        end
    end
end
