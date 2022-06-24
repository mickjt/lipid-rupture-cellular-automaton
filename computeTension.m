% ----------------------------------------------------------------------------
% computeTension: Compute tension in pinned cells
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function cellStatus = computeTension(cellStatus, cellLocation, circleCurrentRadius, numCells, circleMaxRadius, tensionAdhesion, tensionMLV, tensionPin)

for cell = 1:numCells
    if ((cellStatus(cell, 4) == 1)&&(cellStatus(cell, 10) == 0)) % if cell is an unborken pin
        expansionTension = (tensionAdhesion-tensionMLV)*cellLocation(cell, 3)/circleCurrentRadius + tensionMLV; % Irep 3.16
        pinningTension = tensionPin*(circleCurrentRadius - cellLocation(cell, 3))/circleMaxRadius; % Abhay
        cellStatus(cell, 7) = expansionTension + pinningTension;
    end
end
