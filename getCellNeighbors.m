% ----------------------------------------------------------------------------
% getCellNeighbors: Find and store the Moore neighborhood of each cell
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function cellNeighbors = getCellNeighbors(cellKey, automationSize, numCells)


cellNeighbors = zeros(numCells, 8);

% Compute border cell beighbors (except corners)
for borderCell = 2:automationSize-1
    leftRow =cellKey(borderCell,1);
    cellNeighbors(leftRow, 1) = -1;
    cellNeighbors(leftRow, 2) = cellKey(borderCell-1,1);
    cellNeighbors(leftRow, 3) = cellKey(borderCell-1,2);
    cellNeighbors(leftRow, 4) = cellKey(borderCell,2);
    cellNeighbors(leftRow, 5) = cellKey(borderCell+1,2);
    cellNeighbors(leftRow, 6) = cellKey(borderCell+1, 1);
    cellNeighbors(leftRow, 7) = -1;
    cellNeighbors(leftRow, 8) = -1;

    rightRow =cellKey(borderCell,automationSize);
    cellNeighbors(rightRow, 1) = cellKey(borderCell-1, automationSize-1);
    cellNeighbors(rightRow, 2) = cellKey(borderCell-1,automationSize);
    cellNeighbors(rightRow, 3) = -1;
    cellNeighbors(rightRow, 4) = -1;
    cellNeighbors(rightRow, 5) = -1;
    cellNeighbors(rightRow, 6) = cellKey(borderCell+1, automationSize);
    cellNeighbors(rightRow, 7) = cellKey(borderCell+1, automationSize-1);
    cellNeighbors(rightRow, 8) = cellKey(borderCell,automationSize-1);

    topRow = cellKey(1, borderCell);
    cellNeighbors(topRow, 1) = -1;
    cellNeighbors(topRow, 2) = -1;
    cellNeighbors(topRow, 3) = -1;
    cellNeighbors(topRow, 4) = cellKey(1, borderCell+1);
    cellNeighbors(topRow, 5) = cellKey(2, borderCell+1);
    cellNeighbors(topRow, 6) = cellKey(2, borderCell);
    cellNeighbors(topRow, 7) = cellKey(2, borderCell-1);
    cellNeighbors(topRow, 8) = cellKey(1, borderCell-1);

    bottomRow = cellKey(automationSize, borderCell);
    cellNeighbors(bottomRow, 1) = cellKey(automationSize-1, borderCell-1);
    cellNeighbors(bottomRow, 1) = cellKey(automationSize-1, borderCell);
    cellNeighbors(bottomRow, 1) = cellKey(automationSize-1, borderCell+1);
    cellNeighbors(bottomRow, 1) = cellKey(automationSize, borderCell+1);
    cellNeighbors(bottomRow, 1) = -1;
    cellNeighbors(bottomRow, 1) = -1;
    cellNeighbors(bottomRow, 1) = -1;
    cellNeighbors(bottomRow, 1) = cellKey(automationSize, borderCell-1);
end

% Compute neighbors of corner cells
cellNeighbors(cellKey(1,1), 1) = -1;
cellNeighbors(cellKey(1,1), 2) = -1;
cellNeighbors(cellKey(1,1), 3) = -1;
cellNeighbors(cellKey(1,1), 4) = cellKey(1,2);
cellNeighbors(cellKey(1,1), 5) = cellKey(2,2);
cellNeighbors(cellKey(1,1), 6) = cellKey(2,1);
cellNeighbors(cellKey(1,1), 7) = -1;
cellNeighbors(cellKey(1,1), 8) = -1;

cellNeighbors(cellKey(1,automationSize), 1) = -1;
cellNeighbors(cellKey(1,automationSize), 2) = -1;
cellNeighbors(cellKey(1,automationSize), 3) = -1;
cellNeighbors(cellKey(1,automationSize), 4) = -1;
cellNeighbors(cellKey(1,automationSize), 5) = -1;
cellNeighbors(cellKey(1,automationSize), 6) = cellKey(2,automationSize);
cellNeighbors(cellKey(1,automationSize), 7) = cellKey(2, automationSize-1);
cellNeighbors(cellKey(1,automationSize), 8) = cellKey(1, automationSize-1);

cellNeighbors(cellKey(automationSize,1), 1) = -1;
cellNeighbors(cellKey(automationSize,1), 2) = cellKey(automationSize-1,1);
cellNeighbors(cellKey(automationSize,1), 3) = cellKey(automationSize-1,2);
cellNeighbors(cellKey(automationSize,1), 4) = cellKey(automationSize,2);
cellNeighbors(cellKey(automationSize,1), 5) = -1;
cellNeighbors(cellKey(automationSize,1), 6) = -1;
cellNeighbors(cellKey(automationSize,1), 7) = -1;
cellNeighbors(cellKey(automationSize,1), 8) = -1;

cellNeighbors(cellKey(automationSize,automationSize), 1) = cellKey(automationSize-1, automationSize-1);
cellNeighbors(cellKey(automationSize,automationSize), 2) = cellKey(automationSize-1, automationSize);
cellNeighbors(cellKey(automationSize,automationSize), 3) = -1;
cellNeighbors(cellKey(automationSize,automationSize), 4) = -1;
cellNeighbors(cellKey(automationSize,automationSize), 5) = -1;
cellNeighbors(cellKey(automationSize,automationSize), 6) = -1;
cellNeighbors(cellKey(automationSize,automationSize), 7) = -1;
cellNeighbors(cellKey(automationSize,automationSize), 8) = cellKey(automationSize, automationSize-1);

% Compute interior node neighbors
for row =2:automationSize-1
    rowAbove = row-1;
    rowMid = row;
    rowBelow = row+1;

    for col= 2:automationSize-1
        colLeft = col-1;
        colMid = col;
        colRight = col+1;

        cellNeighbors(cellKey(row,col),1) = cellKey(rowAbove,colLeft);
        cellNeighbors(cellKey(row,col),2) = cellKey(rowAbove,colMid);
        cellNeighbors(cellKey(row,col),3) = cellKey(rowAbove,colRight);
        cellNeighbors(cellKey(row,col),4) = cellKey(rowMid,colLeft);
        cellNeighbors(cellKey(row,col),5) = cellKey(rowMid,colRight);
        cellNeighbors(cellKey(row,col),6) = cellKey(rowBelow,colLeft);
        cellNeighbors(cellKey(row,col),7) = cellKey(rowBelow,colMid);
        cellNeighbors(cellKey(row,col),8) = cellKey(rowBelow,colRight);
    end
end
